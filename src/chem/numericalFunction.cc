/*
    File: numericalFunction.cc
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

#include <cando/chem/numericalFunction.h>
#include <math.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>


namespace chem {




#define	ROUND(x) ((x)>0.0)?floor((x)+0.5):-floor((-x)+0.5)



//
// Constructor
//

void	NumericalFunction_O::initialize()
{
    this->Base::initialize();
    this->title = _Nil<core::T_O>();
    this->xAxisName = _Nil<core::T_O>();
    this->yAxisName = _Nil<core::T_O>();
    this->values.clear();
}


//
// copy constructor
//
    NumericalFunction_O::NumericalFunction_O( const NumericalFunction_O& ss ) : core::CxxObject_O(ss)
{
    this->xAxisName = ss.xAxisName;
    this->yAxisName = ss.yAxisName;
    this->xStart = ss.xStart;
    this->xInc = ss.xInc;
    this->values = ss.values;
}

//
// Destructor
//

/*! Return the index of the entry with the largest x-value smaller than or equal to
    the value of (x)
 */
int	NumericalFunction_O::findIndex(double x)
{
int	index;

    if ( this->values.size()== 0 ) {
	return 0;
    }
    if ( this->xInc == 0 ) {
	return 0;
    }
    index = (int)ROUND((double)((x - this->xStart)/this->xInc));
    return index;
}


double	NumericalFunction_O::xValue(int index)
{
double	x;
    x = this->xStart+this->xInc*index;
    return x;
}




void	NumericalFunction_O::setYValueAtIndex(int i, double v)
{
    if ( i < 0 ) {
      for ( int zzz = 0; zzz<(-i); zzz++ ) {
	this->values.insert( this->values.begin(), 0.0 );
      }
      this->xStart = this->xStart + (i*this->xInc);
      this->values[0] = v;
    } else if ( i >= (int)(this->values.size()) ) {
      for ( int zzz = 0; zzz<(i-this->values.size()+1); zzz++ ) {
        this->values.insert( this->values.end(), 0.0 );
      }
      this->values[this->values.size()-1] = v;
    } else {
	this->values[i] = v;
    }
}


double	NumericalFunction_O::getYValueAtIndex(int i)
{
    if ( i < 0 ) {
	return 0.0;
    } else if ( i >= (int)(this->values.size()) ) {
	return 0.0;
    }
    return this->values[i];
}




void	NumericalFunction_O::dumpToXYFile(string fn, string comment)
{
stringstream			dataStream;
FILE*				fOut;
    fOut = fopen(fn.c_str(),"w");
    fprintf( fOut, "# %s %s\n", core::_rep_(this->xAxisName).c_str(), core::_rep_(this->yAxisName).c_str() );
    fprintf( fOut, "# %s\n", comment.c_str() );
    for (auto dp=this->values.begin(); dp!=this->values.end(); dp++ ) {
	fprintf( fOut, "%10.9lf %10.9lf\n", this->xValue(dp-this->values.begin()), *dp );
    }
    fprintf( fOut, "\n");
    fprintf( fOut, "\n");
    fclose(fOut);
}

void	NumericalFunction_O::appendToXYFile(string fn, string comment )
{
stringstream			dataStream;
FILE*				fOut;
    fOut = fopen(fn.c_str(),"a");
    fprintf( fOut, "# %s %s\n", core::_rep_(this->xAxisName).c_str(), core::_rep_(this->yAxisName).c_str() );
    fprintf( fOut, "# %s\n", comment.c_str() );
    for (auto dp=this->values.begin(); dp!=this->values.end(); dp++ ) {
	fprintf( fOut, "%10.9lf %10.9lf\n", this->xValue(dp-this->values.begin()), *dp );
    }
    fprintf( fOut, "\n");
    fprintf( fOut, "\n");
    fclose(fOut);
}








double	NumericalFunction_O::getLowX()
{
    return this->xStart;
}


double	NumericalFunction_O::getHighX()
{
    return this->xStart+this->xInc*this->values.size();
}

double	NumericalFunction_O::getLowY()
{
double	lowY;
    lowY = 99999999.0;
    for (auto dp=this->values.begin(); dp!=this->values.end(); dp++ ) {
	if ( lowY > *dp ) lowY = *dp;
    }
    return lowY;
}

double	NumericalFunction_O::getHighY()
{
double	highY;
    highY = -99999999.0;
    for (auto dp=this->values.begin(); dp!=this->values.end(); dp++ ) {
	if ( highY < *dp ) highY = *dp;
    }
    return highY;
}





/*! Read the text file consisting of lines with v1 v2 v3 v4 v5 ...
    and create a NumericalFunction using x#,y# where # is xindex,yindex
*/
void	NumericalFunction_O::parseFromXYFile(string fn, int xindex, int yindex)
{_OF();
  FILE*				fIn;
  double				x, y;
  double				values[10];
  char				line[256];
  int				n;
  gctools::Vec0<double>			xv, yv;
  this->reset();
  fIn = fopen(fn.c_str(),"r");
  while ( !feof(fIn) ) {
    fgets( line, sizeof(line), fIn );
    n = sscanf( line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &values[0], &values[1], &values[2], &values[3], &values[4],
                &values[5], &values[6], &values[7], &values[8], &values[9] );
    if ( n == 0 ) continue;
    x = 0.0;
    y = 0.0;
    if ( xindex-1 < n ) {
      x = values[xindex-1];
    }  else {
      SIMPLE_ERROR(BF("Index beyond the end of line" ));
    }
    if ( yindex-1 < n ) {
      y = values[yindex-1];
    }  else {
      SIMPLE_ERROR(BF("Index beyond the end of line" ));
    }
    xv.push_back(x);
    yv.push_back(y);
  }
  fclose(fIn);
  this->xInc = xv[1]-xv[0];
  this->xStart = ROUND(xv[0]/this->xInc)*this->xInc;
  this->values.clear();
  for( auto dp=yv.begin(); dp!=yv.end(); dp++ ) {
    this->values.push_back(*dp);
  }
}


void	NumericalFunction_O::reset()
{
    this->values.clear();
}



NumericalFunction_sp NumericalFunction_O::add(NumericalFunction_sp fn)
{_OF();
NumericalFunction_sp		newFunc;
double				xMin, xMax, x;
int				thisIndex, fnIndex, vals;
    if ( this->xInc != fn->xInc ) {
	stringstream ss;
	ss << "xInc values must be identical this->xInc("
		<< this->xInc << ") != fn->xInc("
		<< fn->xInc << ")";
	SIMPLE_ERROR(BF("%s")%ss.str());
    }
    if ( this->values.size() == 0 ) {
	xMin = fn->xStart;
	xMax = fn->xStart+fn->xInc*(fn->values.size());
    } else if ( fn->values.size() == 0 ) {
	xMin = this->xStart;
	xMax = this->xStart+fn->xInc*(fn->values.size());
    } else {
        xMin = fmin(this->xStart,fn->xStart);
        xMax = fmax(this->xStart+this->xInc*(this->values.size()),
		fn->xStart+fn->xInc*(fn->values.size()));
    }
    vals = (int)((xMax-xMin)/this->xInc);
    newFunc = NumericalFunction_O::create();
    newFunc->setXStart(xMin);
    newFunc->setXInc(this->xInc);
    newFunc->values.resize(vals);
    for ( auto ii=newFunc->values.begin(); ii<newFunc->values.end(); ++ii ) {
      *ii = 0.0;
    }
    LOG(BF("resulting NumericalFunction has %d vals") % (newFunc->values.size()) );
    for ( auto dp=newFunc->values.begin(); dp!=newFunc->values.end(); dp++ ) {
	x = newFunc->getXValueAtIndex(dp-newFunc->values.begin());
	LOG(BF("Adding at x=%lf") % (x ) );
	thisIndex = this->findIndex(x);
	LOG(BF("   thisIndex = %d value=%lf") % (thisIndex) % (this->getYValueAtIndex(thisIndex) ) );
	fnIndex = fn->findIndex(x);
	LOG(BF("     fnIndex = %d valut=%lf") % (fnIndex) % (fn->getYValueAtIndex(fnIndex) ) );
	*dp = this->getYValueAtIndex(thisIndex) + fn->getYValueAtIndex(fnIndex);
	LOG(BF("  Wrote value(%lf) to index(%d)") % (*dp) % (dp-newFunc->values.begin()) );
    }
    return newFunc;
}

NumericalFunction_sp NumericalFunction_O::sub(NumericalFunction_sp fn)
{_OF();
NumericalFunction_sp		newFunc;
double				xMin, xMax, x;
int				thisIndex, fnIndex, vals;
    if ( this->xInc != fn->xInc ) {
	stringstream ss;
	ss << "xInc values must be identical this->xInc("
		<< this->xInc << ") != fn->xInc("
		<< fn->xInc << ")";
	SIMPLE_ERROR(BF("%s")%ss.str());
    }
    if ( this->values.size() == 0 ) {
	xMin = fn->xStart;
	xMax = fn->xStart+fn->xInc*(fn->values.size());
    } else if ( fn->values.size() == 0 ) {
	xMin = this->xStart;
	xMax = this->xStart+fn->xInc*(fn->values.size());
    } else {
        xMin = fmin(this->xStart,fn->xStart);
        xMax = fmax(this->xStart+this->xInc*(this->values.size()),
		fn->xStart+fn->xInc*(fn->values.size()));
    }
    vals = (int)((xMax-xMin)/this->xInc);
    newFunc = NumericalFunction_O::create();
    newFunc->setXStart(xMin);
    newFunc->setXInc(this->xInc);
    newFunc->values.resize(vals);
    for (auto ii=newFunc->values.begin(); ii<newFunc->values.end(); ii++ ) *ii = 0.0;
    LOG(BF("resulting NumericalFunction has %d vals") % (newFunc->values.size()) );
    for ( auto dp=newFunc->values.begin(); dp!=newFunc->values.end(); dp++ ) {
	x = newFunc->getXValueAtIndex(dp-newFunc->values.begin());
	LOG(BF("Adding at x=%lf") % (x ) );
	thisIndex = this->findIndex(x);
	LOG(BF("   thisIndex = %d value=%lf") % (thisIndex) % (this->getYValueAtIndex(thisIndex) ) );
	fnIndex = fn->findIndex(x);
	LOG(BF("     fnIndex = %d valut=%lf") % (fnIndex) % (fn->getYValueAtIndex(fnIndex) ) );
	*dp = this->getYValueAtIndex(thisIndex) - fn->getYValueAtIndex(fnIndex);
	LOG(BF("  Wrote value(%lf) to index(%d)") % (*dp) % (dp-newFunc->values.begin()) );
    }
    return newFunc;
}


void	NumericalFunction_O::addValueToHistogram(double x)
{
int		index;
double		cnt, rnd;
    LOG(BF("addValueToHistogram value = %lf    xInc = %lf") % (x) % (this->xInc ) );
    if ( this->values.size() == 0 ) {
	rnd = ROUND(x/this->xInc);
	LOG(BF("addValueToHistogram rnd = %lf") % (rnd ) );
	this->xStart = this->xInc*rnd;
	LOG(BF("addValueToHistogram xStart = %lf") % (this->xStart ) );
	this->values.push_back(1.0);
	return;
    }
    index = this->findIndex(x);
    LOG(BF("addValueToHistogram index= %d") % (index ) );
    cnt = this->getYValueAtIndex(index) + 1.0;
    this->setYValueAtIndex(index,cnt);
}




double	NumericalFunction_O::rmsDifference(NumericalFunction_sp fn)
{
NumericalFunction_sp		newFunc;
double				rms;

    newFunc = this->sub(fn);
    rms = 0.0;
    for ( auto dp=newFunc->values.begin(); dp!=newFunc->values.end(); dp++ ) {
	rms += (*dp)*(*dp);
    }
    rms /= this->values.size();
    return sqrt(rms);
}

void	NumericalFunction_O::scaleXValues(double scale)
{
    this->xStart *= scale;
    this->xInc *= scale;

}

void	NumericalFunction_O::scaleYValues(double scale)
{

    for (auto ip=this->values.begin(); ip!=this->values.end(); ip++ ) {
	*ip *= scale;
    }
}

void	NumericalFunction_O::offsetYValues(double offset)
{
  for (auto ip=this->values.begin(); ip!=this->values.end(); ip++ ) {
    *ip += offset;
  }
}

};
