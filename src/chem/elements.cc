/*
    File: elements.cc
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
#include <clasp/core/symbolToEnumConverter.h>
#include <clasp/core/environment.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/array.h>
#include <clasp/core/ql.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>
#include <cando/geom/color.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/hashTable.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/hashTableEql.h>
#include <cando/chem/chemPackage.h>
#include <cando/chem/elements.h>





namespace chem
{

SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_elementsInfo_PLUS_);
SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_elementToSymbolConverter_PLUS_);
SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_hybridizationToSymbolConverter_PLUS_);
SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_elementColors_PLUS_);
    
SYMBOL_EXPORT_SC_(ChemKwPkg,C);
SYMBOL_EXPORT_SC_(ChemKwPkg,N);
SYMBOL_EXPORT_SC_(ChemKwPkg,O);
SYMBOL_EXPORT_SC_(ChemKwPkg,S);
SYMBOL_EXPORT_SC_(ChemKwPkg,P);     
SYMBOL_EXPORT_SC_(ChemKwPkg,F);
SYMBOL_EXPORT_SC_(ChemKwPkg,Cl);
SYMBOL_EXPORT_SC_(ChemKwPkg,Br);
SYMBOL_EXPORT_SC_(ChemKwPkg,I);


void ElementsInfo_O::initialize() {
  this->_elementFromAtomicSymbol = core::HashTableEq_O::create_default();
};
  
ElementsInfo_sp elementsInfo() {
  ElementsInfo_sp elementsInfo = gc::As<ElementsInfo_sp>(chem::_sym__PLUS_elementsInfo_PLUS_->symbolValue());
  return elementsInfo;
};
  
void set_atomic_info(Element element,const string& as, const string& name, int an, double aw,int am=9999 ) {
  AtomicInfo ai;
  ASSERT(element<element_MAX);
  ai._Valid=true;
  ai._AtomicSymbol = chemkw_intern(as);
  ai._ElementEnum=element;
  ai._ElementName= core::Symbol_O::create_from_string(name);
  ai._AtomicNumber=an;
  ai._IntegerAtomicMass = am;
  ai._RelativeAtomicWeight=aw;
  ai._VdwRadius=1.6 /* Uknown */;
  ai._RealElement = true; 
  ElementsInfo_sp ei = elementsInfo();
  ei->_atomicInfo[element] = ai;
  // Set an element for each atomic number
  if (ei->_atomicNumberToAtomicInfoIndex[an] == -1) {
    ei->_atomicNumberToAtomicInfoIndex[an] = element;
  }
  if (am < ei->_atomicMassToAtomicInfoIndex.size() ) {
    if (ei->_atomicMassToAtomicInfoIndex[am] == -1) {
      ei->_atomicMassToAtomicInfoIndex[am] = element;
    }
  }
};

void    _defineAtomicInfoMapIfNotDefined()
{
  if (chem::_sym__PLUS_elementsInfo_PLUS_->boundP()) {
    printf("%s:%d:%s Returning immediatly\n", __FILE__, __LINE__, __FUNCTION__ );
    return;
  }
  AtomicInfo	invalid;
  invalid._Valid = false;
  auto ei = gctools::GC<ElementsInfo_O>::allocate_with_default_constructor();
  chem::_sym__PLUS_elementsInfo_PLUS_->defparameter(ei);
  if (!chem::_sym__PLUS_elementsInfo_PLUS_.boundp()) {
    printf("%s:%d:%s It's unbound\n", __FILE__, __LINE__, __FUNCTION__ );
    return;
  }
  ei->_atomicInfo.resize(((int)(element_MAX)),invalid);
  ei->_atomicNumberToAtomicInfoIndex.resize((int)element_MAX,-1);
  ei->_atomicMassToAtomicInfoIndex.resize((int)element_mass_MAX,-1);
#define SET_ATOMIC_INFO_FAKE_ELEMENT(element,as,name,an,aw,am)  set_atomic_info(element,as,name,an,aw,am); ei->_atomicInfo[element]._RealElement = false;
#define	SET_VDW_RADIUS(element,vdw) { ei->_atomicInfo[element]._VdwRadius = vdw;};
#define	SET_COVALENT_RADIUS(element,cr) { ei->_atomicInfo[element]._CovalentRadius = cr;};

  LOG(BF("Defining atomicNumberMaps") );
  SET_ATOMIC_INFO_FAKE_ELEMENT(element_LP,"LP","LONE-PAIR",0,0.0,0);
  SET_ATOMIC_INFO_FAKE_ELEMENT(element_Undefined,"UE","UNDEFINED-ELEMENT",0,0.000,0);
  SET_ATOMIC_INFO_FAKE_ELEMENT(element_None,"NE","NO-ELEMENT",0,0.000,0);
  SET_ATOMIC_INFO_FAKE_ELEMENT(element_Wildcard,"WC","WILDCARD-ELEMENT",0,0.00,0);
  SET_ATOMIC_INFO_FAKE_ELEMENT(element_Dummy,"DU","DUMMY-ELEMENT",0,0.000,0);
  set_atomic_info(element_H,"H","Hydrogen",1,1.00794,1);
  set_atomic_info(element_D,"D","Deuterium",1,2.0141017780,2);
  set_atomic_info(element_T,"T","Tritium",1,3.0160492675,3);
  set_atomic_info(element_He,"He","Helium",2,4.002602,4);
  set_atomic_info(element_Li,"Li","Lithium",3,6.941);
  set_atomic_info(element_Be,"Be","Beryllium",4,9.012182);
  set_atomic_info(element_B,"B","Boron",5,10.811,10);
  set_atomic_info(element_C,"C","Carbon",6,12.0107,12);
  set_atomic_info(element_N,"N","Nitrogen",7,14.0067,14);
  set_atomic_info(element_O,"O","Oxygen",8,15.9994,16);
  set_atomic_info(element_F,"F","Fluorine",9,18.9984032);
  set_atomic_info(element_Ne,"Ne","Neon",10,20.1797);
  set_atomic_info(element_Na,"Na","Sodium",11,22.989770);
  set_atomic_info(element_Mg,"Mg","Magnesium",12,24.3050);
  set_atomic_info(element_Al,"Al","Aluminum",13,26.981538);
  set_atomic_info(element_Si,"Si","Silicon",14,28.0855);
  set_atomic_info(element_P,"P","Phosphorus",15,30.973761);
  set_atomic_info(element_S,"S","Sulfur",16,32.065);
  set_atomic_info(element_Cl,"Cl","Chlorine",17,35.453);
  set_atomic_info(element_Ar,"Ar","Argon",18,39.948);
  set_atomic_info(element_K,"K","Potassium",19,39.0983);
  set_atomic_info(element_Ca,"Ca","Calcium",20,40.078);
  set_atomic_info(element_Sc,"Sc","Scandium",21,44.955910);
  set_atomic_info(element_Ti,"Ti","Titanium",22,47.867);
  set_atomic_info(element_V,"V","Vanadium",23,50.9415);
  set_atomic_info(element_Cr,"Cr","Chromium",24,51.9961);
  set_atomic_info(element_Mn,"Mn","Manganese",25,54.938049);
  set_atomic_info(element_Fe,"Fe","Iron",26,55.845);
  set_atomic_info(element_Co,"Co","Cobalt",27,58.933200);
  set_atomic_info(element_Ni,"Ni","Nickel",28,58.6934);
  set_atomic_info(element_Cu,"Cu","Copper",29,63.546);
  set_atomic_info(element_Zn,"Zn","Zinc",30,65.409);
  set_atomic_info(element_Ga,"Ga","Gallium",31,69.723);
  set_atomic_info(element_Ge,"Ge","Germanium",32,72.64);
  set_atomic_info(element_As,"As","Arsenic",33,74.92160);
  set_atomic_info(element_Se,"Se","Selenium",34,78.96);
  set_atomic_info(element_Br,"Br","Bromine",35,79.904);
  set_atomic_info(element_Kr,"Kr","Krypton",36,83.798);
  set_atomic_info(element_Rb,"Rb","Rubidium",37,85.4678);
  set_atomic_info(element_Sr,"Sr","Strontium",38,87.62);
  set_atomic_info(element_Y,"Y","Yttrium",39,88.90585);
  set_atomic_info(element_Zr,"Zr","Zirconium",40,91.224);
  set_atomic_info(element_Nb,"Nb","Niobium",41,92.90638);
  set_atomic_info(element_Mo,"Mo","Molybdenum",42,95.94);
  set_atomic_info(element_Tc,"Tc","Technetium",43,96.906365);
  set_atomic_info(element_Tc,"Tc","Technetium",43,98);
  set_atomic_info(element_Ru,"Ru","Ruthenium",44,101.07);
  set_atomic_info(element_Rh,"Rh","Rhodium",45,102.90550);
  set_atomic_info(element_Pd,"Pd","Palladium",46,106.42);
  set_atomic_info(element_Ag,"Ag","Silver",47,107.8682);
  set_atomic_info(element_Cd,"Cd","Cadmium",48,112.411);
  set_atomic_info(element_In,"In","Indium",49,114.818);
  set_atomic_info(element_Sn,"Sn","Tin",50,118.710);
  set_atomic_info(element_Sb,"Sb","Antimony",51,121.760);
  set_atomic_info(element_Te,"Te","Tellurium",52,127.60);
  set_atomic_info(element_I,"I","Iodine",53,126.90447);
  set_atomic_info(element_Xe,"Xe","Xenon",54,131.293);
  set_atomic_info(element_Cs,"Cs","Cesium",55,132.90545);
  set_atomic_info(element_Ba,"Ba","Barium",56,137.327);
  set_atomic_info(element_La,"La","Lanthanum",57,138.9055);
  set_atomic_info(element_Ce,"Ce","Cerium",58,140.116);
  set_atomic_info(element_Pr,"Pr","Praseodymium",59,140.90765);
  set_atomic_info(element_Nd,"Nd","Neodymium",60,144.24);
  set_atomic_info(element_Sm,"Sm","Samarium",62,150.36);
  set_atomic_info(element_Eu,"Eu","Europium",63,151.964);
  set_atomic_info(element_Gd,"Gd","Gadolinium",64,157.25);
  set_atomic_info(element_Tb,"Tb","Terbium",65,158.92534);
  set_atomic_info(element_Dy,"Dy","Dysprosium",66,162.500);
  set_atomic_info(element_Ho,"Ho","Holmium",67,164.93032);
  set_atomic_info(element_Er,"Er","Erbium",68,167.259);
  set_atomic_info(element_Tm,"Tm","Thulium",69,168.93421);
  set_atomic_info(element_Yb,"Yb","Ytterbium",70,173.04);
  set_atomic_info(element_Lu,"Lu","Lutetium",71,174.967);
  set_atomic_info(element_Hf,"Hf","Hafnium",72,178.49);
  set_atomic_info(element_Ta,"Ta","Tantalum",73,180.9479);
  set_atomic_info(element_W,"W","Wolfram",74,183.84);
  set_atomic_info(element_Re,"Re","Rhenium",75,186.207);
  set_atomic_info(element_Os,"Os","Osmium",76,190.23);
  set_atomic_info(element_Ir,"Ir","Iridium",77,192.217);
  set_atomic_info(element_Pt,"Pt","Platinum",78,195.078);
  set_atomic_info(element_Au,"Au","Gold",79,196.96655);
  set_atomic_info(element_Hg,"Hg","Mercury",80,200.59);
  set_atomic_info(element_Tl,"Tl","Thallium",81,204.3833);
  set_atomic_info(element_Pb,"Pb","Lead",82,207.2);
  set_atomic_info(element_Bi,"Bi","Bismuth",83,208.98038);
  set_atomic_info(element_Th,"Th","Thorium",90,232.0381);
  set_atomic_info(element_Pa,"Pa","Protactinium",91,231.03588);
  set_atomic_info(element_U,"U","Uranium",92,238.02891);
  set_atomic_info(element_Pu,"Pu","Plutonium",94,244.0);
  set_atomic_info(element_Ra,"Ra","Radium",88,226.0);
  SET_VDW_RADIUS(element_H, 1.09);
  SET_VDW_RADIUS(element_D, 1.09);
  SET_VDW_RADIUS(element_T, 1.09);
  SET_VDW_RADIUS(element_C, 1.70);
  SET_VDW_RADIUS(element_N, 1.55);
  SET_VDW_RADIUS(element_O, 1.52);
  SET_VDW_RADIUS(element_F, 1.47);
  SET_VDW_RADIUS(element_P, 1.80);
  SET_VDW_RADIUS(element_S, 1.80);
  SET_VDW_RADIUS(element_Cl, 1.75);
  SET_VDW_RADIUS(element_Cu, 1.4);
  SET_COVALENT_RADIUS(element_H,0.31);
  SET_COVALENT_RADIUS(element_D,0.31);
  SET_COVALENT_RADIUS(element_T,0.31);
  SET_COVALENT_RADIUS(element_He,0.28);
  SET_COVALENT_RADIUS(element_Li,1.28);
  SET_COVALENT_RADIUS(element_Be,0.96);
  SET_COVALENT_RADIUS(element_B,0.84);
  SET_COVALENT_RADIUS(element_C,0.76);
  SET_COVALENT_RADIUS(element_N,0.71);
  SET_COVALENT_RADIUS(element_O,0.66);
  SET_COVALENT_RADIUS(element_F,0.57);
  SET_COVALENT_RADIUS(element_Ne,0.58);
  SET_COVALENT_RADIUS(element_Na,1.66);
  SET_COVALENT_RADIUS(element_Mg,1.41);
  SET_COVALENT_RADIUS(element_Al,1.21);
  SET_COVALENT_RADIUS(element_Si,1.11);
  SET_COVALENT_RADIUS(element_P,1.07);
  SET_COVALENT_RADIUS(element_S,1.05);
  SET_COVALENT_RADIUS(element_Cl,1.02);
  SET_COVALENT_RADIUS(element_Ar,1.06);
  SET_COVALENT_RADIUS(element_K,2.03);
  SET_COVALENT_RADIUS(element_Ca,1.76);
  SET_COVALENT_RADIUS(element_Sc,1.70);
  SET_COVALENT_RADIUS(element_Ti,1.60);
  SET_COVALENT_RADIUS(element_V,1.53);
  SET_COVALENT_RADIUS(element_Cr,1.39);
  SET_COVALENT_RADIUS(element_Mn,1.61);
  SET_COVALENT_RADIUS(element_Fe,1.52);
  SET_COVALENT_RADIUS(element_Co,1.50);
  SET_COVALENT_RADIUS(element_Ni,1.24);
  SET_COVALENT_RADIUS(element_Cu,1.32);
  SET_COVALENT_RADIUS(element_Zn,1.22);
  SET_COVALENT_RADIUS(element_Ga,1.22);
  SET_COVALENT_RADIUS(element_Ge,1.20);
  SET_COVALENT_RADIUS(element_As,1.19);
  SET_COVALENT_RADIUS(element_Se,1.20);
  SET_COVALENT_RADIUS(element_Br,1.20);
  SET_COVALENT_RADIUS(element_Kr,1.16);
  SET_COVALENT_RADIUS(element_Rb,2.20);
  SET_COVALENT_RADIUS(element_Sr,1.95);
  SET_COVALENT_RADIUS(element_Y,1.90);
  SET_COVALENT_RADIUS(element_Zr,1.75);
  SET_COVALENT_RADIUS(element_Nb,1.64);
  SET_COVALENT_RADIUS(element_Mo,1.54);
  SET_COVALENT_RADIUS(element_Tc,1.47);
  SET_COVALENT_RADIUS(element_Ru,1.46);
  SET_COVALENT_RADIUS(element_Rh,1.42);
  SET_COVALENT_RADIUS(element_Pd,1.39);
  SET_COVALENT_RADIUS(element_Ag,1.45);
  SET_COVALENT_RADIUS(element_Cd,1.44);
  SET_COVALENT_RADIUS(element_In,1.42);
  SET_COVALENT_RADIUS(element_Sn,1.39);
  SET_COVALENT_RADIUS(element_Sb,1.39);
  SET_COVALENT_RADIUS(element_Te,1.38);
  SET_COVALENT_RADIUS(element_I,1.39);
  SET_COVALENT_RADIUS(element_Xe,1.40);
  SET_COVALENT_RADIUS(element_Cs,2.44);
  SET_COVALENT_RADIUS(element_Ba,2.15);
  SET_COVALENT_RADIUS(element_La,2.07);
  SET_COVALENT_RADIUS(element_Ce,2.04);
  SET_COVALENT_RADIUS(element_Pr,2.03);
  SET_COVALENT_RADIUS(element_Nd,2.01);
  SET_COVALENT_RADIUS(element_Pm,1.99);
  SET_COVALENT_RADIUS(element_Sm,1.98);
  SET_COVALENT_RADIUS(element_Eu,1.98);
  SET_COVALENT_RADIUS(element_Gd,1.96);
  SET_COVALENT_RADIUS(element_Tb,1.94);
  SET_COVALENT_RADIUS(element_Dy,1.92);
  SET_COVALENT_RADIUS(element_Ho,1.92);
  SET_COVALENT_RADIUS(element_Er,1.89);
  SET_COVALENT_RADIUS(element_Tm,1.90);
  SET_COVALENT_RADIUS(element_Yb,1.87);
  SET_COVALENT_RADIUS(element_Lu,1.87);
  SET_COVALENT_RADIUS(element_Hf,1.75);
  SET_COVALENT_RADIUS(element_Ta,1.70);
  SET_COVALENT_RADIUS(element_W,1.62);
  SET_COVALENT_RADIUS(element_Re,1.51);
  SET_COVALENT_RADIUS(element_Os,1.44);
  SET_COVALENT_RADIUS(element_Ir,1.41);
  SET_COVALENT_RADIUS(element_Pt,1.36);
  SET_COVALENT_RADIUS(element_Au,1.36);
  SET_COVALENT_RADIUS(element_Hg,1.32);
  SET_COVALENT_RADIUS(element_Tl,1.45);
  SET_COVALENT_RADIUS(element_Pb,1.46);
  SET_COVALENT_RADIUS(element_Bi,1.48);
  SET_COVALENT_RADIUS(element_Po,1.40);
  SET_COVALENT_RADIUS(element_At,1.50);
  SET_COVALENT_RADIUS(element_Rn,1.50);
  SET_COVALENT_RADIUS(element_Fr,2.60);
  SET_COVALENT_RADIUS(element_Ra,2.21);
  SET_COVALENT_RADIUS(element_Ac,2.15);
  SET_COVALENT_RADIUS(element_Th,2.06);
  SET_COVALENT_RADIUS(element_Pa,2.00);
  SET_COVALENT_RADIUS(element_U,1.96);
  SET_COVALENT_RADIUS(element_Pu,1.87);

  ei->_hybridizationInfo.resize(((int)(hybridization_MAX)));
  ei->_hybridizationInfo[hybridization_undefined]._HybridizationEnum = hybridization_undefined;
  ei->_hybridizationInfo[hybridization_undefined]._HybridizationSymbol = chemkw_intern("SP?");
  ei->_hybridizationInfo[hybridization_sp]._HybridizationEnum = hybridization_sp;
  ei->_hybridizationInfo[hybridization_sp]._HybridizationSymbol = chemkw_intern("SP");
  ei->_hybridizationInfo[hybridization_sp2]._HybridizationEnum = hybridization_sp2;
  ei->_hybridizationInfo[hybridization_sp2]._HybridizationSymbol = chemkw_intern("SP2");
  ei->_hybridizationInfo[hybridization_sp3]._HybridizationEnum = hybridization_sp3;
  ei->_hybridizationInfo[hybridization_sp3]._HybridizationSymbol = chemkw_intern("SP3");
  for ( auto ai=ei->_atomicInfo.begin(); ai!=ei->_atomicInfo.end(); ai++ )
  {
    if ( ai->_Valid ) {
      ei->_elementFromAtomicSymbol->setf_gethash(ai->_AtomicSymbol, core::make_fixnum(ai->_ElementEnum));
    }
  }
}


void initializeElementsAndHybridization()
{
  _defineAtomicInfoMapIfNotDefined();
  core::SymbolToEnumConverter_sp elementConverter = core::SymbolToEnumConverter_O::create("Element");
  ElementsInfo_sp ei = elementsInfo();
  for ( auto ai=ei->_atomicInfo.begin(); ai!=ei->_atomicInfo.end(); ai++ )
  {
    if ( ai->_Valid )
    {
      core::Symbol_sp sym = ai->_AtomicSymbol;
      elementConverter->addSymbolEnumPair(sym,ai->_AtomicSymbol,ai->_ElementEnum);
    }
  }
  _sym__PLUS_elementToSymbolConverter_PLUS_->defconstant(elementConverter);

    //
    // Now generate the hybridization symbols
    //
  core::SymbolToEnumConverter_sp hybridizationConverter = core::SymbolToEnumConverter_O::create("Hybridization");
  for ( auto hi=ei->_hybridizationInfo.begin(); hi!=ei->_hybridizationInfo.end(); hi++ )
  {
    core::Symbol_sp sym = hi->_HybridizationSymbol;
    hybridizationConverter->addSymbolEnumPair(sym,hi->_HybridizationSymbol,hi->_HybridizationEnum);
  }
  _sym__PLUS_hybridizationToSymbolConverter_PLUS_->defconstant(hybridizationConverter);

//  core::af_def(ChemPkg,"elementFromAtomNameString",&elementFromAtomNameString);
  core::Symbol_sp elementColors = _sym__PLUS_elementColors_PLUS_;
  core::HashTableEql_sp dict = core::HashTableEql_O::create_default();
  dict->setf_gethash(core::clasp_make_fixnum(element_C),geom::Color_O::systemColor(kw::_sym_cyan));
  dict->setf_gethash(core::clasp_make_fixnum(element_H),geom::Color_O::systemColor(kw::_sym_white));
  dict->setf_gethash(core::clasp_make_fixnum(element_O),geom::Color_O::systemColor(kw::_sym_red));
  dict->setf_gethash(core::clasp_make_fixnum(element_N),geom::Color_O::systemColor(kw::_sym_blue));
  dict->setf_gethash(core::clasp_make_fixnum(element_S),geom::Color_O::systemColor(kw::_sym_yellow));
  dict->setf_gethash(core::clasp_make_fixnum(element_P),geom::Color_O::systemColor(kw::_sym_yellow));
  dict->setf_gethash(core::clasp_make_fixnum(element_Undefined),geom::Color_O::systemColor(kw::_sym_magenta));
  elementColors->defparameter(dict);
}



/*!
 * Return the symbol associated with the element
 */
core::Symbol_sp symbolFromElement(Element element)
{
  core::SymbolToEnumConverter_sp elementConverter = _sym__PLUS_elementToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
  return elementConverter->symbolForEnum<Element>(element);
}

Element elementForSymbol(core::Symbol_sp sym)
{
  core::SymbolToEnumConverter_sp elementConverter(_sym__PLUS_elementToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>());
  return elementConverter->enumForSymbol<chem::Element>(sym);
}



bool elementIsRealElement(Element element)
{
  HARD_ASSERT(element>0 && element<element_MAX);
  ElementsInfo_sp ei = elementsInfo();
  AtomicInfo& ai = ei->_atomicInfo[element];
  HARD_ASSERT(ai._Valid);
  return ai._RealElement;
}


Hybridization hybridizationForSymbol(core::Symbol_sp sym)
{
  core::SymbolToEnumConverter_sp hybridizationConverter = _sym__PLUS_hybridizationToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
  return hybridizationConverter->enumForSymbol<Hybridization>(sym);
}

/*!
 * Return the symbol associated with the hybridization
 */
DOCGROUP(cando)
CL_DEFUN core::Symbol_sp symbolFromHybridization(chem::Hybridization hybridization)
{
  core::SymbolToEnumConverter_sp hybridizationConverter = _sym__PLUS_hybridizationToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
  return hybridizationConverter->symbolForEnum<Hybridization>(hybridization);
}


CL_DOCSTRING(R"dx(Return hybridizaton value for string)dx")
DOCGROUP(cando)
CL_DEFUN chem::Hybridization hybridizationForString(const string& sym)
{
  core::SymbolToEnumConverter_sp hybridizationConverter = _sym__PLUS_hybridizationToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
  return hybridizationConverter->enumForSymbol<Hybridization>(_lisp->intern(sym));
}


DOCGROUP(cando)
CL_DEFUN core::Symbol_sp hybridizationSymbolFromHybridization(chem::Hybridization h)
{
  core::SymbolToEnumConverter_sp hybridizationConverter = _sym__PLUS_hybridizationToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
  return hybridizationConverter->symbolForEnum<Hybridization>(h);
}

DOCGROUP(cando)
CL_DEFUN core::Symbol_sp atomicSymbolFromElement(chem::Element h)
{
  core::SymbolToEnumConverter_sp elementConverter = _sym__PLUS_elementToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
  return elementConverter->symbolForEnum<Element>(h);
}





uint maxTotalBondOrderForElement(Element element)
{
  if ( element == element_C ) return 4;
  if ( element == element_N ) return 4;
  if ( element == element_O ) return 3;
  if ( element == element_H ) return 1;
  if ( element == element_S ) return 6;
  if ( element == element_P ) return 6;
  return 4;
}

CL_DOCSTRING(R"dx(Return VdW radius for element)dx")
DOCGROUP(cando)
CL_DEFUN double vdwRadiusForElement(::chem::Element element)
{
  ElementsInfo_sp ei = elementsInfo();
  return ei->_atomicInfo[element]._VdwRadius;
}

bool elementIsHydrogen(Element e)
{
  if ( e == element_H ) return true;
  if ( e == element_D ) return true;
  if ( e == element_T ) return true;
  return false;
}

Element elementFromSymbol(core::Symbol_sp sym)
{
  core::SymbolToEnumConverter_sp elementConverter = _sym__PLUS_elementToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
  Element element = elementConverter->enumForSymbol<Element>(sym);
  return element;
}


double maximumBondLengthBetweenElements(Element ie1, Element ie2 )
{
  ElementsInfo_sp ei = elementsInfo();
  double cr1 = ei->_atomicInfo[ie1]._CovalentRadius;
  double cr2 = ei->_atomicInfo[ie2]._CovalentRadius;
  double maxLength = cr1+cr2+0.56;
  return maxLength;
}

char convertCase(char c, bool caseInsensitive)
{
  if (caseInsensitive) {
    return toupper(c);
  }
  return c;
}

string isolateElementName(const string& name, bool caseInsensitive)
{
  char firstChar = convertCase(name[0],caseInsensitive);
  if (name == "" ) {
    SIMPLE_ERROR(BF("Cannot isolate element name from empty string"));
  } else if (name.size() == 1) {
    goto ONE_CHAR_ELEMENT;
  }
  // Names like C_xxx C.xxx etc the first char is the element
  if (!isalpha(name[1])) {
    goto ONE_CHAR_ELEMENT;
  }
  // If name is like Ca.xxx Nb.yyy the first two characters are element
  if (name.size()>2&&!isalnum(name[2])) {
    string element = "  ";
    element[0] = firstChar;
    element[1] = name[1];
    return element;
  }
  // If first char upper and second char is lower case name as in CaXXX
  if (isalpha(name[0])&&isalpha(name[1])&&isupper(name[0])&&islower(name[1])) {
    string element = "  ";
    element[0] = firstChar;
    element[1] = name[1];
    return element;
  }
  // If first two characters are uppercase as in CAXXX
  // then use first character only
  if (isupper(name[0])&&isupper(name[1])) {
    goto ONE_CHAR_ELEMENT;
  }
  // If first two characters are lowercase as in caxxx
  // then use first character only
  if (isupper(name[0])&&isupper(name[1])) {
    goto ONE_CHAR_ELEMENT;
  }
  core::eval::funcall(cl::_sym_warn,
                      core::Str_O::create("Cannot identify element from ~a"),
                      core::Str_O::create(name));
 ONE_CHAR_ELEMENT:
  string element = " ";
  element[0] = firstChar;
  return element;
}
                      

  
Element elementFromAtomNameStringBasic(const std::string& name, bool caseInsensitive)
{
  ElementsInfo_sp ei = elementsInfo();
  string elementName = isolateElementName(name,caseInsensitive);
  core::Symbol_sp se = chemkw_intern(elementName);
  core::T_mv val = ei->_elementFromAtomicSymbol->gethash(se);

  if (val.notnilp()) {
    Element element = (Element)val.unsafe_fixnum();
    if (ei->_atomicInfo[element]._RealElement) {
      return element;
    }
  }

  core::eval::funcall(cl::_sym_warn,
                      core::Str_O::create("Could not determine element from name ~a case-insensitive(~a)"),
                      core::Str_O::create(name),
                      _lisp->_boolean(caseInsensitive));
  return element_None;
}



DOCGROUP(cando)
CL_DEFUN ::chem::Element elementFromAtomNameStringCaseInsensitive(const std::string& name)
{
    return elementFromAtomNameStringBasic(name,true);
}

DOCGROUP(cando)
CL_DEFUN ::chem::Element elementFromAtomNameString(const string& name)
{
    return elementFromAtomNameStringBasic(name,false);
}

DOCGROUP(cando)
CL_DEFUN core::List_sp all_element_symbols() {
  ElementsInfo_sp ei = elementsInfo();
  ql::list ll;
  for ( auto ai=ei->_atomicInfo.begin(); ai!=ei->_atomicInfo.end(); ai++ ) {
    if ( ai->_Valid ) {
      core::Symbol_sp sym = ai->_AtomicSymbol;
      ll << sym;
    }
  }
  return ll.cons();
}

  
DOCGROUP(cando)
CL_DEFUN chem::Element elementForAtomicNumber(int atomicNumber)
{
  ElementsInfo_sp ei = elementsInfo();
  if (atomicNumber<0 || atomicNumber>element_MAX || ei->_atomicNumberToAtomicInfoIndex[atomicNumber] ==-1 ) {
    SIMPLE_ERROR(BF("Could not determine element for atomic number %d") % atomicNumber);
  }
  return ei->_atomicInfo[ei->_atomicNumberToAtomicInfoIndex[atomicNumber]]._ElementEnum;
}

DOCGROUP(cando)
CL_DEFUN chem::Element elementForIntegerAtomicMass(int atomicMass)
{
  ElementsInfo_sp ei = elementsInfo();
  if (atomicMass<0 || atomicMass>element_mass_MAX || ei->_atomicMassToAtomicInfoIndex[atomicMass] ==-1 ) {
    SIMPLE_ERROR(BF("Could not determine element for integer atomic mass %d") % atomicMass);
  }
  return ei->_atomicInfo[ei->_atomicMassToAtomicInfoIndex[atomicMass]]._ElementEnum;
}

DOCGROUP(cando)
CL_DEFUN uint atomicNumberForElement(chem::Element element)
{
  ElementsInfo_sp ei = elementsInfo();
  HARD_ASSERT(element>0 && element<element_MAX);
  AtomicInfo& ai = ei->_atomicInfo[element];
  HARD_ASSERT(ai._Valid);
  return ai._AtomicNumber;
}

DOCGROUP(cando)
CL_DEFUN uint integerAtomicMassForElement(chem::Element element)
{
  ElementsInfo_sp ei = elementsInfo();
  HARD_ASSERT(element>0 && element<element_MAX);
  AtomicInfo& ai = ei->_atomicInfo[element];
  HARD_ASSERT(ai._Valid);
  return ai._IntegerAtomicMass;
}

DOCGROUP(cando)
CL_DEFUN double atomicWeightForElement(chem::Element element)
{
  ElementsInfo_sp ei = elementsInfo();
    HARD_ASSERT(element>0 && element<element_MAX);
    AtomicInfo& ai = ei->_atomicInfo[element];
    HARD_ASSERT(ai._Valid);
    return ai._RelativeAtomicWeight;
}




};
