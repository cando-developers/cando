#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/symbolToEnumConverter.h"
#include "core/environment.h"
#include "core/binder.h"
#include "core/symbolTable.h"
#include "core/wrappers.h"
#include "addon/color.h"
#include "addon/symbolTable.h"
#include "symbolTable.h"
#include "addon/integerKeyObjectDictionary.h"
#include "chemPackage.h"
#include "elements.h"





namespace chem
{

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

struct	AtomicInfo {
    bool	_Valid;
    bool	_RealElement;
    Element	_ElementEnum;
    string	_AtomicSymbol;
    string	_ElementName;
    int		_AtomicNumber;
    double	_RelativeAtomicWeight;
    double	_VdwRadius;
    double	_CovalentRadius;
};


struct HybridizationInfo
{
    Hybridization	_HybridizationEnum;
    string		_HybridizationSymbol;
};

bool    atomicInfoDefined = false;
vector<AtomicInfo> 		atomicInfo;
map<string,Element>		elementFromAtomicSymbol;
set<string>			twoCharacterElementNames;
set<string>			oneCharacterElementNames;
vector<HybridizationInfo>	hybridizationInfo;


void    _defineAtomicInfoMapIfNotDefined()
{_G();
    if ( atomicInfoDefined ) return;
    AtomicInfo	invalid;
    invalid._Valid = false;
    atomicInfo.clear();
    elementFromAtomicSymbol.clear();
    atomicInfo.resize(((int)(element_MAX)),invalid);

#define SET_ATOMIC_INFO(element,as,name,an,aw) {AtomicInfo ai;\
					ASSERT(element<element_MAX);\
					ai._Valid=true;\
    					ai._AtomicSymbol=as;\
					ai._ElementEnum=element;\
					ai._ElementName=name;\
					ai._AtomicNumber=an;\
					ai._RelativeAtomicWeight=aw;\
					ai._VdwRadius=1.0 /* Uknown */;\
					ai._RealElement = true; \
				       	atomicInfo[element] = ai;\
				    };
#define SET_ATOMIC_INFO_FAKE_ELEMENT(element,as,name,an,aw)  SET_ATOMIC_INFO(element,as,name,an,aw); atomicInfo[element]._RealElement = false;
#define	SET_VDW_RADIUS(element,vdw) { atomicInfo[element]._VdwRadius = vdw;};
#define	SET_COVALENT_RADIUS(element,cr) { atomicInfo[element]._CovalentRadius = cr;};

    atomicInfoDefined = true;
    LOG(BF("Defining atomicNumberMaps") );
    SET_ATOMIC_INFO_FAKE_ELEMENT(element_LP,"LP","LONE-PAIR",0,0.0);
    SET_ATOMIC_INFO_FAKE_ELEMENT(element_Undefined,"UE","UNDEFINED-ELEMENT",0,0.0000);
    SET_ATOMIC_INFO_FAKE_ELEMENT(element_None,"NE","NO-ELEMENT",0,0.0000);
    SET_ATOMIC_INFO_FAKE_ELEMENT(element_Wildcard,"WC","WILDCARD-ELEMENT",0,0.0000);
    SET_ATOMIC_INFO_FAKE_ELEMENT(element_Dummy,"DU","DUMMY-ELEMENT",0,0.0000);
    SET_ATOMIC_INFO(element_H,"H","Hydrogen",1,1.00794);
    SET_ATOMIC_INFO(element_D,"D","Deuterium",1,2.0141017780);
    SET_ATOMIC_INFO(element_T,"T","Tritium",1,3.0160492675);
    SET_ATOMIC_INFO(element_He,"He","Helium",2,4.002602);
    SET_ATOMIC_INFO(element_Li,"Li","Lithium",3,6.941);
    SET_ATOMIC_INFO(element_Be,"Be","Beryllium",4,9.012182);
    SET_ATOMIC_INFO(element_B,"B","Boron",5,10.811);
    SET_ATOMIC_INFO(element_C,"C","Carbon",6,12.0107);
    SET_ATOMIC_INFO(element_N,"N","Nitrogen",7,14.0067);
    SET_ATOMIC_INFO(element_O,"O","Oxygen",8,15.9994);
    SET_ATOMIC_INFO(element_F,"F","Fluorine",9,18.9984032);
    SET_ATOMIC_INFO(element_Ne,"Ne","Neon",10,20.1797);
    SET_ATOMIC_INFO(element_Na,"Na","Sodium",11,22.989770);
    SET_ATOMIC_INFO(element_Mg,"Mg","Magnesium",12,24.3050);
    SET_ATOMIC_INFO(element_Al,"Al","Aluminum",13,26.981538);
    SET_ATOMIC_INFO(element_Si,"Si","Silicon",14,28.0855);
    SET_ATOMIC_INFO(element_P,"P","Phosphorus",15,30.973761);
    SET_ATOMIC_INFO(element_S,"S","Sulfur",16,32.065);
    SET_ATOMIC_INFO(element_Cl,"Cl","Chlorine",17,35.453);
    SET_ATOMIC_INFO(element_Ar,"Ar","Argon",18,39.948);
    SET_ATOMIC_INFO(element_K,"K","Potassium",19,39.0983);
    SET_ATOMIC_INFO(element_Ca,"Ca","Calcium",20,40.078);
    SET_ATOMIC_INFO(element_Sc,"Sc","Scandium",21,44.955910);
    SET_ATOMIC_INFO(element_Ti,"Ti","Titanium",22,47.867);
    SET_ATOMIC_INFO(element_V,"V","Vanadium",23,50.9415);
    SET_ATOMIC_INFO(element_Cr,"Cr","Chromium",24,51.9961);
    SET_ATOMIC_INFO(element_Mn,"Mn","Manganese",25,54.938049);
    SET_ATOMIC_INFO(element_Fe,"Fe","Iron",26,55.845);
    SET_ATOMIC_INFO(element_Co,"Co","Cobalt",27,58.933200);
    SET_ATOMIC_INFO(element_Ni,"Ni","Nickel",28,58.6934);
    SET_ATOMIC_INFO(element_Cu,"Cu","Copper",29,63.546);
    SET_ATOMIC_INFO(element_Zn,"Zn","Zinc",30,65.409);
    SET_ATOMIC_INFO(element_Ga,"Ga","Gallium",31,69.723);
    SET_ATOMIC_INFO(element_Ge,"Ge","Germanium",32,72.64);
    SET_ATOMIC_INFO(element_As,"As","Arsenic",33,74.92160);
    SET_ATOMIC_INFO(element_Se,"Se","Selenium",34,78.96);
    SET_ATOMIC_INFO(element_Br,"Br","Bromine",35,79.904);
    SET_ATOMIC_INFO(element_Kr,"Kr","Krypton",36,83.798);
    SET_ATOMIC_INFO(element_Rb,"Rb","Rubidium",37,85.4678);
    SET_ATOMIC_INFO(element_Sr,"Sr","Strontium",38,87.62);
    SET_ATOMIC_INFO(element_Y,"Y","Yttrium",39,88.90585);
    SET_ATOMIC_INFO(element_Zr,"Zr","Zirconium",40,91.224);
    SET_ATOMIC_INFO(element_Nb,"Nb","Niobium",41,92.90638);
    SET_ATOMIC_INFO(element_Mo,"Mo","Molybdenum",42,95.94);
    SET_ATOMIC_INFO(element_Tc,"Tc","Technetium",43,96.906365);
    SET_ATOMIC_INFO(element_Tc,"Tc","Technetium",43,98);
    SET_ATOMIC_INFO(element_Ru,"Ru","Ruthenium",44,101.07);
    SET_ATOMIC_INFO(element_Rh,"Rh","Rhodium",45,102.90550);
    SET_ATOMIC_INFO(element_Pd,"Pd","Palladium",46,106.42);
    SET_ATOMIC_INFO(element_Ag,"Ag","Silver",47,107.8682);
    SET_ATOMIC_INFO(element_Cd,"Cd","Cadmium",48,112.411);
    SET_ATOMIC_INFO(element_In,"In","Indium",49,114.818);
    SET_ATOMIC_INFO(element_Sn,"Sn","Tin",50,118.710);
    SET_ATOMIC_INFO(element_Sb,"Sb","Antimony",51,121.760);
    SET_ATOMIC_INFO(element_Te,"Te","Tellurium",52,127.60);
    SET_ATOMIC_INFO(element_I,"I","Iodine",53,126.90447);
    SET_ATOMIC_INFO(element_Xe,"Xe","Xenon",54,131.293);
    SET_ATOMIC_INFO(element_Cs,"Cs","Cesium",55,132.90545);
    SET_ATOMIC_INFO(element_Ba,"Ba","Barium",56,137.327);
    SET_ATOMIC_INFO(element_La,"La","Lanthanum",57,138.9055);
    SET_ATOMIC_INFO(element_Ce,"Ce","Cerium",58,140.116);
    SET_ATOMIC_INFO(element_Pr,"Pr","Praseodymium",59,140.90765);
    SET_ATOMIC_INFO(element_Nd,"Nd","Neodymium",60,144.24);
    SET_ATOMIC_INFO(element_Sm,"Sm","Samarium",62,150.36);
    SET_ATOMIC_INFO(element_Eu,"Eu","Europium",63,151.964);
    SET_ATOMIC_INFO(element_Gd,"Gd","Gadolinium",64,157.25);
    SET_ATOMIC_INFO(element_Tb,"Tb","Terbium",65,158.92534);
    SET_ATOMIC_INFO(element_Dy,"Dy","Dysprosium",66,162.500);
    SET_ATOMIC_INFO(element_Ho,"Ho","Holmium",67,164.93032);
    SET_ATOMIC_INFO(element_Er,"Er","Erbium",68,167.259);
    SET_ATOMIC_INFO(element_Tm,"Tm","Thulium",69,168.93421);
    SET_ATOMIC_INFO(element_Yb,"Yb","Ytterbium",70,173.04);
    SET_ATOMIC_INFO(element_Lu,"Lu","Lutetium",71,174.967);
    SET_ATOMIC_INFO(element_Hf,"Hf","Hafnium",72,178.49);
    SET_ATOMIC_INFO(element_Ta,"Ta","Tantalum",73,180.9479);
    SET_ATOMIC_INFO(element_W,"W","Wolfram",74,183.84);
    SET_ATOMIC_INFO(element_Re,"Re","Rhenium",75,186.207);
    SET_ATOMIC_INFO(element_Os,"Os","Osmium",76,190.23);
    SET_ATOMIC_INFO(element_Ir,"Ir","Iridium",77,192.217);
    SET_ATOMIC_INFO(element_Pt,"Pt","Platinum",78,195.078);
    SET_ATOMIC_INFO(element_Au,"Au","Gold",79,196.96655);
    SET_ATOMIC_INFO(element_Hg,"Hg","Mercury",80,200.59);
    SET_ATOMIC_INFO(element_Tl,"Tl","Thallium",81,204.3833);
    SET_ATOMIC_INFO(element_Pb,"Pb","Lead",82,207.2);
    SET_ATOMIC_INFO(element_Bi,"Bi","Bismuth",83,208.98038);
    SET_ATOMIC_INFO(element_Th,"Th","Thorium",90,232.0381);
    SET_ATOMIC_INFO(element_Pa,"Pa","Protactinium",91,231.03588);
    SET_ATOMIC_INFO(element_U,"U","Uranium",92,238.02891);
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
SET_COVALENT_RADIUS(element_H,0.31)
SET_COVALENT_RADIUS(element_D,0.31)
SET_COVALENT_RADIUS(element_T,0.31)
SET_COVALENT_RADIUS(element_He,0.28)
SET_COVALENT_RADIUS(element_Li,1.28)
SET_COVALENT_RADIUS(element_Be,0.96)
SET_COVALENT_RADIUS(element_B,0.84)
SET_COVALENT_RADIUS(element_C,0.76)
SET_COVALENT_RADIUS(element_N,0.71)
SET_COVALENT_RADIUS(element_O,0.66)
SET_COVALENT_RADIUS(element_F,0.57)
SET_COVALENT_RADIUS(element_Ne,0.58)
SET_COVALENT_RADIUS(element_Na,1.66)
SET_COVALENT_RADIUS(element_Mg,1.41)
SET_COVALENT_RADIUS(element_Al,1.21)
SET_COVALENT_RADIUS(element_Si,1.11)
SET_COVALENT_RADIUS(element_P,1.07)
SET_COVALENT_RADIUS(element_S,1.05)
SET_COVALENT_RADIUS(element_Cl,1.02)
SET_COVALENT_RADIUS(element_Ar,1.06)
SET_COVALENT_RADIUS(element_K,2.03)
SET_COVALENT_RADIUS(element_Ca,1.76)
SET_COVALENT_RADIUS(element_Sc,1.70)
SET_COVALENT_RADIUS(element_Ti,1.60)
SET_COVALENT_RADIUS(element_V,1.53)
SET_COVALENT_RADIUS(element_Cr,1.39)
SET_COVALENT_RADIUS(element_Mn,1.61)
SET_COVALENT_RADIUS(element_Fe,1.52)
SET_COVALENT_RADIUS(element_Co,1.50)
SET_COVALENT_RADIUS(element_Ni,1.24)
SET_COVALENT_RADIUS(element_Cu,1.32)
SET_COVALENT_RADIUS(element_Zn,1.22)
SET_COVALENT_RADIUS(element_Ga,1.22)
SET_COVALENT_RADIUS(element_Ge,1.20)
SET_COVALENT_RADIUS(element_As,1.19)
SET_COVALENT_RADIUS(element_Se,1.20)
SET_COVALENT_RADIUS(element_Br,1.20)
SET_COVALENT_RADIUS(element_Kr,1.16)
SET_COVALENT_RADIUS(element_Rb,2.20)
SET_COVALENT_RADIUS(element_Sr,1.95)
SET_COVALENT_RADIUS(element_Y,1.90)
SET_COVALENT_RADIUS(element_Zr,1.75)
SET_COVALENT_RADIUS(element_Nb,1.64)
SET_COVALENT_RADIUS(element_Mo,1.54)
SET_COVALENT_RADIUS(element_Tc,1.47)
SET_COVALENT_RADIUS(element_Ru,1.46)
SET_COVALENT_RADIUS(element_Rh,1.42)
SET_COVALENT_RADIUS(element_Pd,1.39)
SET_COVALENT_RADIUS(element_Ag,1.45)
SET_COVALENT_RADIUS(element_Cd,1.44)
SET_COVALENT_RADIUS(element_In,1.42)
SET_COVALENT_RADIUS(element_Sn,1.39)
SET_COVALENT_RADIUS(element_Sb,1.39)
SET_COVALENT_RADIUS(element_Te,1.38)
SET_COVALENT_RADIUS(element_I,1.39)
SET_COVALENT_RADIUS(element_Xe,1.40)
SET_COVALENT_RADIUS(element_Cs,2.44)
SET_COVALENT_RADIUS(element_Ba,2.15)
SET_COVALENT_RADIUS(element_La,2.07)
SET_COVALENT_RADIUS(element_Ce,2.04)
SET_COVALENT_RADIUS(element_Pr,2.03)
SET_COVALENT_RADIUS(element_Nd,2.01)
SET_COVALENT_RADIUS(element_Pm,1.99)
SET_COVALENT_RADIUS(element_Sm,1.98)
SET_COVALENT_RADIUS(element_Eu,1.98)
SET_COVALENT_RADIUS(element_Gd,1.96)
SET_COVALENT_RADIUS(element_Tb,1.94)
SET_COVALENT_RADIUS(element_Dy,1.92)
SET_COVALENT_RADIUS(element_Ho,1.92)
SET_COVALENT_RADIUS(element_Er,1.89)
SET_COVALENT_RADIUS(element_Tm,1.90)
SET_COVALENT_RADIUS(element_Yb,1.87)
SET_COVALENT_RADIUS(element_Lu,1.87)
SET_COVALENT_RADIUS(element_Hf,1.75)
SET_COVALENT_RADIUS(element_Ta,1.70)
SET_COVALENT_RADIUS(element_W,1.62)
SET_COVALENT_RADIUS(element_Re,1.51)
SET_COVALENT_RADIUS(element_Os,1.44)
SET_COVALENT_RADIUS(element_Ir,1.41)
SET_COVALENT_RADIUS(element_Pt,1.36)
SET_COVALENT_RADIUS(element_Au,1.36)
SET_COVALENT_RADIUS(element_Hg,1.32)
SET_COVALENT_RADIUS(element_Tl,1.45)
SET_COVALENT_RADIUS(element_Pb,1.46)
SET_COVALENT_RADIUS(element_Bi,1.48)
SET_COVALENT_RADIUS(element_Po,1.40)
SET_COVALENT_RADIUS(element_At,1.50)
SET_COVALENT_RADIUS(element_Rn,1.50)
SET_COVALENT_RADIUS(element_Fr,2.60)
SET_COVALENT_RADIUS(element_Ra,2.21)
SET_COVALENT_RADIUS(element_Ac,2.15)
SET_COVALENT_RADIUS(element_Th,2.06)
SET_COVALENT_RADIUS(element_Pa,2.00)
SET_COVALENT_RADIUS(element_U,1.96)

    hybridizationInfo.resize(((int)(hybridization_MAX)));
    hybridizationInfo[hybridization_undefined]._HybridizationEnum = hybridization_undefined;
    hybridizationInfo[hybridization_undefined]._HybridizationSymbol = "sp?";
    hybridizationInfo[hybridization_sp]._HybridizationEnum = hybridization_sp;
    hybridizationInfo[hybridization_sp]._HybridizationSymbol = "sp";
    hybridizationInfo[hybridization_sp2]._HybridizationEnum = hybridization_sp2;
    hybridizationInfo[hybridization_sp2]._HybridizationSymbol = "sp2";
    hybridizationInfo[hybridization_sp3]._HybridizationEnum = hybridization_sp3;
    hybridizationInfo[hybridization_sp3]._HybridizationSymbol = "sp3";
    for ( vector<AtomicInfo>::iterator ai=atomicInfo.begin(); ai!=atomicInfo.end(); ai++ )
    {
	if ( ai->_Valid )
	{
	    elementFromAtomicSymbol[ai->_AtomicSymbol] = ai->_ElementEnum;
	    if ( ai->_Valid && ai->_RealElement )
	    {
		if ( ai->_AtomicSymbol.size() == 2 )
		{
		    twoCharacterElementNames.insert(ai->_AtomicSymbol);
		} else if ( ai->_AtomicSymbol.size() == 1 )
		{
		    oneCharacterElementNames.insert(ai->_AtomicSymbol);
		}
	    }
	}
    }
}


void initializeElementsAndHybridization()
{_G();
    _defineAtomicInfoMapIfNotDefined();
    core::SymbolToEnumConverter_sp elementConverter = core::SymbolToEnumConverter_O::create("Element");
    for ( vector<AtomicInfo>::iterator ai=atomicInfo.begin(); ai!=atomicInfo.end(); ai++ )
    {
	if ( ai->_Valid )
	{
	    core::Symbol_sp sym = chemkw_intern(ai->_AtomicSymbol);
	    FIX_ME(); elementConverter->addSymbolEnumPair(sym,chemkw_intern(ai->_AtomicSymbol),ai->_ElementEnum);
	}
    }
    _sym__PLUS_elementToSymbolConverter_PLUS_->defconstant(elementConverter);

    //
    // Now generate the hybridization symbols
    //
    core::SymbolToEnumConverter_sp hybridizationConverter = core::SymbolToEnumConverter_O::create("Hybridization");
    for ( vector<HybridizationInfo>::iterator hi=hybridizationInfo.begin(); hi!=hybridizationInfo.end(); hi++ )
    {
	core::Symbol_sp sym = chemkw_intern(hi->_HybridizationSymbol);
	hybridizationConverter->addSymbolEnumPair(sym,_lisp->intern(hi->_HybridizationSymbol),hi->_HybridizationEnum);
    }
    _sym__PLUS_hybridizationToSymbolConverter_PLUS_->defconstant(hybridizationConverter);

    core::af_def(ChemPkg,"elementFromAtomNameString",&elementFromAtomNameString);
    core::af_def(ChemPkg,"elementFromAtomNameStringCaseInsensitive",&elementFromAtomNameStringCaseInsensitive);


    core::Symbol_sp elementColors = _sym__PLUS_elementColors_PLUS_;
    addon::IntegerKeyObjectDictionary_sp dict = addon::IntegerKeyObjectDictionary_O::create();
    _lisp->defvar(elementColors,dict);
    dict->set(element_C,addon::Color_O::systemColor(kw::_sym_cyan));
    dict->set(element_H,addon::Color_O::systemColor(kw::_sym_white));
    dict->set(element_O,addon::Color_O::systemColor(kw::_sym_red));
    dict->set(element_N,addon::Color_O::systemColor(kw::_sym_blue));
    dict->set(element_S,addon::Color_O::systemColor(kw::_sym_yellow));
    dict->set(element_P,addon::Color_O::systemColor(kw::_sym_yellow));
    dict->set(element_Undefined,addon::Color_O::systemColor(kw::_sym_magenta));
}



/*!
 * Return the symbol associated with the element
 */
core::Symbol_sp symbolFromElement(Element element)
{_G();
    core::SymbolToEnumConverter_sp elementConverter = _sym__PLUS_elementToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
    return elementConverter->symbolForEnum<Element>(element);
}

Element elementForSymbol(core::Symbol_sp sym)
{_G();
    core::SymbolToEnumConverter_sp elementConverter(_sym__PLUS_elementToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>());
    return elementConverter->enumForSymbol<chem::Element>(sym);
}



bool elementIsRealElement(Element element)
{
    HARD_ASSERT(element>0 && element<element_MAX);
    AtomicInfo& ai = atomicInfo[element];
    HARD_ASSERT(ai._Valid);
    return ai._RealElement;
}


Hybridization hybridizationForSymbol(core::Symbol_sp sym)
{_G();
    core::SymbolToEnumConverter_sp hybridizationConverter = _sym__PLUS_hybridizationToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
    return hybridizationConverter->enumForSymbol<Hybridization>(sym);
}

/*!
 * Return the symbol associated with the hybridization
 */
core::Symbol_sp symbolFromHybridization(Hybridization hybridization)
{_G();
    core::SymbolToEnumConverter_sp hybridizationConverter = _sym__PLUS_hybridizationToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
    return hybridizationConverter->symbolForEnum<Hybridization>(hybridization);
}



Hybridization hybridizationForString(const string& sym)
{_G();
    core::SymbolToEnumConverter_sp hybridizationConverter = _sym__PLUS_hybridizationToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
    return hybridizationConverter->enumForSymbol<Hybridization>(_lisp->intern(sym));
}


core::Symbol_sp hybridizationSymbolFromHybridization(Hybridization h)
{_G();
    core::SymbolToEnumConverter_sp hybridizationConverter = _sym__PLUS_hybridizationToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
    return hybridizationConverter->symbolForEnum<Hybridization>(h);
}

core::Symbol_sp	atomicSymbolFromElement(Element h)
{_G();
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
    return 4;
}

double vdwRadiusForElement(Element element)
{
    return atomicInfo[element]._VdwRadius;
}

    bool elementIsHydrogen(Element e)
    {
	if ( e == element_H ) return true;
	if ( e == element_D ) return true;
	if ( e == element_T ) return true;
	return false;
    }

Element elementFromSymbol(core::Symbol_sp sym)
{_G();
    core::SymbolToEnumConverter_sp elementConverter = _sym__PLUS_elementToSymbolConverter_PLUS_->symbolValue().as<core::SymbolToEnumConverter_O>();
    Element element = elementConverter->enumForSymbol<Element>(sym);
    return element;
}


    double maximumBondLengthBetweenElements(Element ie1, Element ie2 )
    {
	double cr1 = atomicInfo[ie1]._CovalentRadius;
	double cr2 = atomicInfo[ie2]._CovalentRadius;
	double maxLength = cr1+cr2+0.56;
	return maxLength;
    }


Element elementFromAtomNameStringBasic(const string& name, bool caseInsensitive)
{_G();
    LOG(BF("elementFromName name(%s)") % name.c_str()  );
    stringstream namePart;
    if ( name.size() == 2 )
    {
	if ( caseInsensitive )
	{
	    namePart << (char)(toupper(name[0])) << (char)(tolower(name[1]));
	} else
	{
	    namePart << name[0] << name[1];
	}
	LOG(BF("Trying two character names against(%s)") % namePart.str() );
	if ( twoCharacterElementNames.count(namePart.str()) > 0 )
	{
	    Element element = elementFromAtomicSymbol[namePart.str()];
	    LOG(BF("Found match element(%d)") % element );
	    return element;
	}
	LOG(BF("Didn't match two character names"));
	return element_None;
    }
    namePart.str("");
    if ( caseInsensitive )
    {
	namePart << (char)(toupper(name[0]));
    } else
    {
	namePart << name[0];
    }
    LOG(BF("Trying one character names against(%s)") % namePart.str() );
    if ( oneCharacterElementNames.count(namePart.str()) )
    {
	Element element = elementFromAtomicSymbol[namePart.str()];
	LOG(BF("Found match element(%d)") % element );
	return element;
    }
    LOG(BF("No match, returning with element=(?)"));
    return element_None;
}



Element elementFromAtomNameStringCaseInsensitive(const string& name)
{
    return elementFromAtomNameStringBasic(name,true);
}

Element elementFromAtomNameString(const string& name)
{
    return elementFromAtomNameStringBasic(name,false);
}


uint atomicNumberForElement(Element element)
{
    HARD_ASSERT(element>0 && element<element_MAX);
    AtomicInfo& ai = atomicInfo[element];
    HARD_ASSERT(ai._Valid);
    return ai._AtomicNumber;
}

double atomicWeightForElement(Element element)
{
    HARD_ASSERT(element>0 && element<element_MAX);
    AtomicInfo& ai = atomicInfo[element];
    HARD_ASSERT(ai._Valid);
    return ai._RelativeAtomicWeight;
}




};
