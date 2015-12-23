       
       




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
//	void	archive(core::ArchiveP node);
private:
	string		title;
	string		xAxisName;
	string		yAxisName;
	double		xInc;
	double		xStart;
    std::vector<double>	values;

public:
	static NumericalFunction_sp create(const string& xaxis,
				const string& yaxis,
				double xstart, double xinc) ;
public:


	void	empty() { this->values.clear();};
	void	setTitle(string x) {this->title = x;};
	void	setXAxisName(string x) {this->xAxisName = x;};
	void	setYAxisName(string y) {this->yAxisName = y;};

	string	getTitle() {return this->title;};
	string	getXAxisName() {return this->xAxisName;};
	string	getYAxisName() {return this->yAxisName;};

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


inline	NumericalFunction_sp NumericalFunction_O::create(const string& xaxis,
				const string& yaxis,
				double xstart, double xinc) 
{
    GC_ALLOCATE(NumericalFunction_O, f ); // RP_Create<NumericalFunction_O>(e);
    f->setXAxisName(xaxis);
    f->setYAxisName(yaxis);
    f->setXStart(xstart);
    f->setXInc(xinc);
    return f;
}

};
TRANSLATE(chem::NumericalFunction_O);
#endif
