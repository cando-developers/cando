/*
    File: nVector.cc
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

//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/foundation.h>
#include <clasp/core/array.h>
#include <cando/chem/nVector.h>

namespace chem {
/*
 *      dotProduct
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the DOT product of two vectors.
 *      return:  (this) . x
 */
double	dotProduct( NVector_sp x, NVector_sp y )
{
  double          dDot;
  ASSERT(x->length()==y->length());
  dDot = 0.0;
  double* px = &(*x)[0];
  double* py = &(*y)[0];
  for ( size_t i(0), iEnd(x->length()); i<iEnd; ++i ) dDot += px[i]*py[i];
  return dDot;
}


double	squared(NVector_sp x)
{
  double          dDot;
  dDot = 0.0;
  double* dp = &(*x)[0];
  for ( size_t i(0), iEnd(x->length()); i<iEnd; ++i ) {
    dDot += dp[i]*dp[i];
  }
  return(dDot);
}

double magnitude(NVector_sp me)
{
  double  dDot;
  dDot = squared(me);
  return(sqrt(dDot));
}


double	rmsMagnitude(NVector_sp me)
{
  ASSERT(me->length()>0);
  double  dDot;
  dDot = squared(me);
  dDot /= (double)me->length();
  return(sqrt(dDot));
}


double	angleWithVector(NVector_sp me, NVector_sp other)
{
#define	VERY_SMALL 1.0e-6
  double	lenThis, lenOther, dot;
  lenThis = magnitude(me);
  lenOther = magnitude(other);
  dot = dotProduct(me->asSmartPtr(),other);
  if ( fabs(lenThis) > VERY_SMALL && fabs(lenOther)>VERY_SMALL) {
    dot /= (lenThis*lenOther);
    if ( dot > 1.0 ) dot = 1.0;
    if ( dot < -1.0 ) dot = -1.0;
    return acos(dot);
  }
  return 0.0;
}

void XPlusYTimesScalar( NVector_sp o, NVector_sp x, NVector_sp y, double s )
{
 LOG(BF("this->_Values.size() = %d") % o->length() );
 LOG(BF("x->size() = %d") % x->length() );
 LOG(BF("y->size() = %d") % y->length() );
 ASSERTP( o->length() == x->length(),"NVector_O::addTimesScalar>>mismatched length");
 ASSERT( o->length() == y->length());
 double* po = &(*o)[0];
 double* px = &(*x)[0];
 double* py = &(*y)[0];
 for ( int i(0),iEnd(o->length()); i<iEnd; ++i ) {
   po[i] = px[i]+py[i]*s;
 }
}

/*
 *      copyVector
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Copy the contents of nvVector into another.
 */
void	copyVector(NVector_sp dest, NVector_sp orig)
{
  ASSERTP(dest->length()==orig->length(),"NVector_O::copy>>mismatched length");
  double* ddest = &(*dest)[0];
  double* dorig = &(*orig)[0];
  for ( size_t i(0),iEnd(orig->length()); i<iEnd; ++i ) {
    ddest[i] = dorig[i];
  }
}

void inPlaceAddTimesScalar( NVector_sp result, NVector_sp dir, double s )
{_OF();
  ASSERT( result->length() == dir->length());
  double* presult = &(*result)[0];
  double* pdir = &(*dir)[0];
  for ( size_t i(0), iEnd(result->length()); i<iEnd; ++i ) {
    presult[i] += pdir[i]*s;
  }
}


double	rmsDistanceFrom(NVector_sp u, NVector_sp v)
{_OF();
  ASSERT(u->length() == v->length());
  double	e, sum;
  sum = 0.0;
  double* du = &(*u)[0];
  double* dv = &(*v)[0];
  for ( size_t i(0), iEnd(u->length()); i<iEnd; ++i ) {
    e = du[i]-dv[i];
    sum += e*e;
  }
  sum /= (float)(u->length());
  return sqrt(sum);
}


};


#if 0

/*
 *      File:   nVector.h
 *
 *      Description:
 *              Maintain N dimensional vectors.
 */
#include <iostream>
#include <math.h>
#include <cando/chem/nVector.h>
#include <clasp/core/wrappers.h>



namespace chem {


/*
 *      NVector_O
 *
 *	Author:	Christian Schafmeister (1991)
 *
 */
    NVector_O::NVector_O( const NVector_O& orig ) : NVector_O::Base(orig)
{
    this->_Size = orig._Size;
    this->_Values.assign(orig._Values.begin(),orig._Values.end());
}



void	NVector_O::initialize()
{
    this->setup(0);
}

void	NVector_O::setup( uint iElements )
{
    this->_Size = iElements;
    this->_Values.resize(iElements,0.0);
    LOG(BF("Created NVector_O@%0X with %d entries") % &(this->_Values) % this->_Size );
}


/*
 *      copy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Copy the contents of nvVector into another.
 */
void	NVector_O::copyIntoMe(NVector_sp orig)
{_OF();
    ASSERTP(this->_Size == orig->_Size,"NVector_O::copy>>mismatched size");
    this->_Values.assign(orig->_Values.begin(),orig->_Values.end());
}


core::T_sp copy_nvector(NVector_sp orig ) const
{ 
  GC_COPY(NVector_O, vv, orig );
  return vv;
}

/*
 *      zero
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the vector with zeros
 */
void	NVector_O::zero()
{
    this->_Values.assign(this->_Size,0.0);
}



/*
 *      fill
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the vector with a number.
 */
CL_LISPIFY_NAME(nvector-fill);
CL_DEFMETHOD void NVector_O::fill( double dValue )
{
    this->_Values.assign(this->_Size,dValue);
}

/*
 *      Return the largest absolute value in this vector.
 */
double	NVector_O::maxAbsValue( )
{
double		maxAbsVal;
vector<double>::iterator	dp;

    maxAbsVal = 0.0;
    for ( dp=this->_Values.begin(); dp!=this->_Values.end(); dp++ ) {
	maxAbsVal= MAX(maxAbsVal,fabs(*dp));
    }
    return maxAbsVal;
}

/*
 *      Return the smallest absolute value in this vector.
 */
double	NVector_O::minAbsValue( )
{
double		minAbsVal;
vector<double>::iterator	dp;
    minAbsVal = 9.0e99;
    for ( dp=this->_Values.begin(); dp!=this->_Values.end(); dp++ ) {
	minAbsVal= MIN(minAbsVal,fabs(*dp));
    }
    return minAbsVal;
}


/*
 *      Return the largest value in this vector.
 */
double	NVector_O::maxValue( )
{
double		maxVal;
vector<double>::iterator	dp;

    maxVal = -9.0e99;
    for ( dp=this->_Values.begin(); dp!=this->_Values.end(); dp++ ) {
	maxVal= MAX(maxVal,*dp);
    }
    return maxVal;
}

/*
 *      Return the smallest absolute value in this vector.
 */
double	NVector_O::minValue( )
{
double		minVal;
vector<double>::iterator	dp;
    minVal = 9.0e99;
    for ( dp=this->_Values.begin(); dp!=this->_Values.end(); dp++ ) {
	minVal= MIN(minVal,*dp);
    }
    return minVal;
}


/*
 *      describe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe the NVector_sp to cout.
 */
void	NVector_O::describe()
{
    uint	i;
    std::cout << "NVector_O length= " << this->_Size << std::endl;
    for ( i=0; i<this->_Size; i++ )
        std::cout << "---[" << i << "] = " << this->_Values[i] << std::endl;
}



string	NVector_O::description() const
{ 
stringstream	ss;
    ss << "NVector[";
    for (uint i=0; i<this->_Size; i++ ) {
	ss << this->_Values[i] << " ";
    }
    ss << "]";
    return ss.str();
}


void	NVector_O::debugDump(const string& name)
{_OF();
    IMPLEMENT_ME();
}


/*
 *      add
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add two vectors and put the result in (this)
 */
CL_LISPIFY_NAME(nvector-dest-x-_PLUS_-y);
CL_DEFMETHOD void NVector_O::dest_x__PLUS__y( NVector_sp x, NVector_sp y )
{_OF();
vector<double>::iterator	resi, xi,yi;
    ASSERTP( this->_Values.size() == x->size(),"NVector_O::add>>mismatched size");
    ASSERTP( this->_Values.size() == y->size(),"NVector_O::add>>mismatched size");
    for ( resi = this->_Values.begin(), xi=x->_Values.begin(), yi=y->_Values.begin();
	    resi!=this->_Values.end(); resi++, xi++, yi++ ) {
	*resi = *xi+*yi;
    }
}


/*!
 *	this = x + y*s
 */

/*
 *      addScalar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a scalar to each value (this)
 */
void	NVector_O::addScalar( NVector_sp x, double s )
{
vector<double>::iterator	resi, xi;
    LOG(BF("this->_Values.size() = %d") % this->_Values.size() );
    LOG(BF("x->size() = %d") % x->size() );
    ASSERTP( this->_Values.size() == x->size(),"NVector_O::addTimesScalar>>mismatched size");
    for ( resi = this->_Values.begin(), xi=x->_Values.begin();
	    resi!=this->_Values.end(); resi++, xi++ ) {
	*resi = *xi+s;
    }
}





/*
 *      sub
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add two vectors and return the result in (this)
 */
void	NVector_O::sub( NVector_sp x, NVector_sp y )
{_OF();
vector<double>::iterator	resi, xi,yi;
    ASSERT( this->_Values.size() == x->size());
    ASSERT( this->_Values.size() == y->size());
    for ( resi = this->_Values.begin(), xi=x->_Values.begin(), yi=y->_Values.begin();
	    resi!=this->_Values.end(); resi++, xi++, yi++ ) {
	*resi = *xi-*yi;
    }
}



/*
 *      timesScalar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Multiply a vector by a scalar and put the result in (this)
 *      A = X * s
 */
void	NVector_O::timesScalar( NVector_sp x, double d )
{
vector<double>::iterator	resi, xi;
    for ( resi=this->_Values.begin(),xi=x->_Values.begin(); resi!=this->_Values.end(); resi++,xi++ ) {
	*resi = (*xi)*d;
    }
}







  double	NVector_O::squared()
  {
    double          dDot;
    vector<double>::iterator	resi;
    dDot = 0.0;
    for ( resi=this->_Values.begin(); resi!=this->_Values.end(); resi++ ) {
      dDot += (*resi)*(*resi);
    }
    return(dDot);
  }



  void	NVector_O::inPlaceTimesScalar( double s)
  {
    vector<double>::iterator	resi;
    for ( resi=this->_Values.begin(); resi!=this->_Values.end(); resi++ ) {
      *resi = (*resi)*s;
    }
  }

#ifdef XML_ARCHIVE
  void	NVector_O::archive(core::ArchiveP node)
  {
    node->archiveVectorDouble( "vec", this->_Values );
    LOG(BF("size of this->_Values=%d") % this->_Values.size()  );
    if ( node->loading() )
    {
      this->_Size = this->_Values.size();
    }
  }
#endif


  void	NVector_O::writeMathematica(const string& fileName)
  {
    std::ofstream	fout;
    uint		y;
    fout.open(fileName.c_str());
    fout << "{" << std::endl;
    fout << std::setiosflags(std::ios::fixed)
         << std::setw(16)
         << std::setprecision(8);
    for ( y = 0; y<this->_Size; y++ ) {
      if ( y!= 0 ) {
        fout << "," << std::endl;
      }
      fout << this->element(y);
    }
    fout << "}" << std::endl;
    fout.close();
  }





};
#endif
