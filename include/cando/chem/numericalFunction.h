/*
    File: numericalFunction.h
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
       
       




#ifndef	NumericalFunction_H
#define NumericalFunction_H

#include <clasp/core/common.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <set>


#include <cando/chem/chemPackage.h>


namespace chem {

inline string	NumericalFunction_XmlName() {return "NumericalFunctionMold";};


typedef	struct	{
	double	x;
	double	y;
} XYPair;


SMART(NumericalFunction);


/*! A NumericalFunction is an X/Y function where the X-range is defined by the (xStart) value
 *and the (xInc) interval length.
 */
SMART(NumericalFunction);
class NumericalFunction_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,NumericalFunction_O,"NumericalFunction",core::CxxObject_O);

public:
	void initialize();
public:
  core::T_sp		title;
  core::T_sp		xAxisName;
  core::T_sp		yAxisName;
  double		xInc;
  double		xStart;
  gctools::Vec0<double>	values;

public:
  static NumericalFunction_sp create(core::T_sp xaxis,
                                     core::T_sp yaxis,
                                     double xstart, double xinc) ;
public:


	void	empty() { this->values.clear();};
  void	setTitle(core::T_sp x) {this->title = x;};
  void	setXAxisName(core::T_sp x) {this->xAxisName = x;};
  void	setYAxisName(core::T_sp y) {this->yAxisName = y;};

  core::T_sp	getTitle() {return this->title;};
  core::T_sp	getXAxisName() {return this->xAxisName;};
  core::T_sp	getYAxisName() {return this->yAxisName;};

	void	setXStart(double xs) {this->xStart = xs;};
	void	setXInc(double inc) {this->xInc = inc; };
	double	getXStart() {return this->xStart;};
	double	getXInc() {return this->xInc; };
	void		reset();

	int	findIndex(double x);
	double	xValue(int i);

	int	getNumberOfValues() { return this->values.size(); };

	double		getLowX();
	double		getHighX();
	double		getLowY();
	double		getHighY();

	double	getYValueAtIndex(int i);
	void	setYValueAtIndex(int i, double d);
	double	getXValueAtIndex(int i) {return this->xStart+this->xInc*i;};
	double	getYValueAtX(double x);
	void	setYValueAtX(double x);

	void	addValueToHistogram(double x);

	void	appendValue(double y) {this->values.push_back(y);};

//	void		addValues(double ax, double ay);
//	double		interpolate(double ax);

	double		rmsDifference(NumericalFunction_sp fn);

	NumericalFunction_sp	add(NumericalFunction_sp fn);
	NumericalFunction_sp	sub(NumericalFunction_sp fn);

	void		scaleXValues(double scale );
	void		scaleYValues(double scale );
	void		offsetYValues(double offset );



	void		parseFromXYFile(string fn,int xoffset,int yoffset);

	void		dumpToXYFile(string xyfn,string comment);
	void		appendToXYFile(string xyfn,string comment);
//    std::ostream&	dumpToStream(ostream& out);


	NumericalFunction_O( const NumericalFunction_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(NumericalFunction_O);
};



// Constructor function


inline	NumericalFunction_sp NumericalFunction_O::create(core::T_sp xaxis,
                                                         core::T_sp yaxis,
                                                         double xstart, double xinc) 
{
  auto  f  = gctools::GC<NumericalFunction_O>::allocate_with_default_constructor();
    f->setXAxisName(xaxis);
    f->setYAxisName(yaxis);
    f->setXStart(xstart);
    f->setXInc(xinc);
    return f;
}

};
#endif
