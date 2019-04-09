/*
    File: elements.h
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
  CL_DEFUN core::Symbol_sp symbolFromElement(chem::Element element);
  CL_DEFUN bool elementIsHydrogen(chem::Element element);

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

  Element elementForAtomicNumber(int atomicNumber);

Element elementForIntegerAtomicMass(int atomicMass); // this is used in SMARTS
    
  uint	atomicNumberForElement(Element element);
  uint	integerAtomicMassForElement(Element element);

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
      if (o.fixnump()) {
        this->_v = chem::elementForAtomicNumber(o.unsafe_fixnum());
        return;
      } else if (core::Symbol_sp sym = o.asOrNull<core::Symbol_O>() ) {
        this->_v = chem::elementForSymbol(sym);
        return;
      }
      SIMPLE_ERROR(BF("Could not convert %s to an element - either provide an element symbol or an atomic number")% _rep_(o));
    };

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
