/*
    File: representedEntityNameSet.h
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
#ifndef	RepresentedEntityNameSet_H
#define RepresentedEntityNameSet_H





#include <cando/chem/entityNameSet.h>


#include <cando/chem/chemPackage.h>


namespace chem {

SMART(MapOfMonomerNamesToAtomIndexers);

SMART(RepresentedEntityNameSet);
class RepresentedEntityNameSet_O : public EntityNameSet_O
{
    LISP_CLASS(chem,ChemPkg,RepresentedEntityNameSet_O,"RepresentedEntityNameSet",EntityNameSet_O);
 public:
    static RepresentedEntityNameSet_sp make(core::Symbol_sp representativeName);
public:
  void	initialize();
private:
	core::Symbol_sp	_Representative;
public:
    virtual RepresentativeList_sp expandedRepresentativeList() const;
public:
    string __repr__() const { return this->description(); };
	string description() const;
	bool	hasRepresentative();
	void	setRepresentative(core::Symbol_sp cc );
	core::Symbol_sp	getRepresentative();

	/*! Return a copy of this object */
	EntityNameSetBase_sp copy() const;


		/*! Return a RepresentedEntityNameSet containing itself as its
		 * contents - this is used when creating trainers
		 * so that the contents of the RepresentedEntityNameSet are not locked
		 * in at the time the trainer is defined.
		 */
//	RepresentedEntityNameSet_sp asUnexpandedRepresentedEntityNameSet();

CL_LISPIFY_NAME("testRepresentedEntityNameSet");
CL_DEFMETHOD 	string	testRepresentedEntityNameSet() { return "testRepresentedEntityNameSet"; };
//	virtual	string	getGeneralKey();

	RepresentedEntityNameSet_O( const RepresentedEntityNameSet_O& ss ); //!< Copy constructor



	DEFAULT_CTOR_DTOR(RepresentedEntityNameSet_O);
};






};


TRANSLATE(chem::RepresentedEntityNameSet_O);
#endif
