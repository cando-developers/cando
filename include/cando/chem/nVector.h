/*
    File: nVector.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *      File:   nVector.h
 *
 *      Description:
 *              Routines to maintain N dimensional vectors
 */

#ifndef NVECTOR_H
#define NVECTOR_H

#include <clasp/core/common.h>
//#include "core/archiveNode.h"



#include <cando/chem/chemPackage.h>


namespace chem {

    SMART(NVector );
    class NVector_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,NVector_O,"NVector",core::CxxObject_O);

public:
	void initialize();
public:
	static NVector_sp create(uint sz);
public:
	uint     	_Size;
	vector<double>	_Values;
        NVector_O( const NVector_O& orig );
	void	setup(uint sz);
	void	describe();
	string	description() const;
	void	debugDump(const string& name );

	void	copyIntoMe(NVector_sp orig);
	core::T_sp	deepCopy() const;

	void	zero();		// NVectorZero();
	void	fill(double d);
	void	dest_x__PLUS__y( NVector_sp nvX, NVector_sp nvY);		// NVectorAdd(a,x,y)
		//! this = x + d*sigma
	void	XPlusYTimesScalar( NVector_sp nvX, NVector_sp nvD, double sigma);

	void	addScalar( NVector_sp nvX, double x);		// NVectorAdd(a,x,y)
		//!  this += nvY*x;
	void	inPlaceAddTimesScalar( NVector_sp nvY, double x);	// NVectorAdd(a,x,y)
	void	sub( NVector_sp nvX, NVector_sp nvY );	// NVectorSub(a,x,y)
	void	timesScalar( NVector_sp nvX, double dS);
	double	dotProduct( NVector_sp nvY);
	double	magnitude();
	double	rmsDistanceFrom(NVector_sp v);
	double	rmsMagnitude();
	double	squared();

		//! Return the angle wrt other
	double	angleWithVector(NVector_sp other);
	void		writeMathematica(const string& name);
//	adapt::QDomNode_sp	termDifferencesAsXml(const string& thisName, NVector_sp other, const string& otherName );
	double	maxAbsValue();
	double	minAbsValue();
	double	maxValue();
	double	minValue();
	uint	size()	{ return this->_Size; };
	double	element(uint i)	{return this->_Values[i];};
	double	getElement(uint i)	{_G(); ASSERT_lt(i,this->_Values.size());return this->_Values[i];};
	void	setElement(uint i, double d) { this->_Values[i]=d;};
	void	addToElement(uint i, double d) { this->_Values[i]+=d;};
//	void	setValues(int s,double*d) {this->iSize=s;this->dPValues=d;};

	void	inPlaceTimesScalar(double s);

//	void	archive(core::ArchiveP node);
	DEFAULT_CTOR_DTOR(NVector_O);
    };


#if 0
inline	NVector_sp	new_NVector_sp() {
    return NVector_sp(new NVector_O(0));
}

inline	NVector_sp	new_NVector_sp(uint size) {
    return NVector_sp(new NVector_O(size));
}

inline	NVector_sp	create_NVector(uint size) {
    return NVector_sp(new NVector_O(size));
}

inline	NVector_sp	new_NVector_sp(const NVector_O& orig) {
    NVector_sp copy = RP_Copy<NVector_O>(&orig);
    return copy;
}

#endif


};

TRANSLATE(chem::NVector_O);
#endif          /* ifndef NVECTOR_H */
