#ifndef	elements_H
#define elements_H

#include <clasp/core/common.h>
#include <clasp/core/symbol.h>

#include <cando/chem/elements.fwd.h>


#include <cando/chem/chemPackage.h>


namespace chem
{



void	initializeElementsAndHybridization();

/*! Return true if the element is a real one and not a dummy or undefined etc
 */
bool elementIsRealElement(Element element);


uint maxTotalBondOrderForElement(Element element);

    Element elementForSymbol(core::Symbol_sp sym);
core::Symbol_sp symbolFromElement(Element element);
bool elementIsHydrogen(Element element);

core::Symbol_sp symbolFromHybridization(Hybridization element);




core::Symbol_sp keywordSymbolForHybridization(Hybridization hyb);

Hybridization hybridizationForSymbol(core::Symbol_sp sym);

double vdwRadiusForElement(Element element);

Element elementFromAtomNameStringBasic(const string& name,bool caseSensitive);

Element elementFromAtomNameStringCaseInsensitive(const string& name);

/*! Case sensitive */
    Element elementFromAtomNameString(const string& name);

	/*! Get the atomic symbol from the Element
	 */
    core::Symbol_sp atomicSymbolFromElement(Element element);

	/*! Get the hybridization symbol string "sp", "sp2", "sp3" 
	 * from the hybridization
	 */
    core::Symbol_sp hybridizationSymbolFromHybridization(Hybridization hyb);


uint	atomicNumberForElement(Element element);

double	atomicWeightForElement(Element element);


double maximumBondLengthBetweenElements(Element element1, Element element2 );




};






namespace translate
{
    template <>
    struct translate::from_object<chem::Element>
    {
	typedef	chem::Element ExpectedType;
	typedef	chem::Element DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
	    core::Symbol_sp sym = o.as<core::Symbol_O>();
	    this->_v = chem::elementForSymbol(sym);
	}
    };



    template <>
    struct	to_object<chem::Element>
    {
	typedef	core::Symbol_sp	ExpectedType;
	typedef	core::Symbol_sp	DeclareType;
	static core::T_sp convert( chem::Element v)
	{_G();
	    return (chem::symbolFromElement(v));
	}
    };






    template <>
    struct translate::from_object<chem::Hybridization>
    {
	typedef	chem::Hybridization ExpectedType;
	typedef	chem::Hybridization DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
	    core::Symbol_sp sym = o.as<core::Symbol_O>();
	    this->_v = chem::hybridizationForSymbol(sym);
	}
    };



    template <>
    struct	to_object<chem::Hybridization>
    {
	typedef	core::Symbol_sp	ExpectedType;
	typedef	core::Symbol_sp	DeclareType;
	static core::T_sp convert(chem::Hybridization v)
	{_G();
	    return (chem::symbolFromHybridization(v));
	}
    };


};




#endif // elements_H
