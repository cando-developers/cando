       
#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//


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

    REGISTER_CLASS(chem,NVector_O);

NVector_sp NVector_O::create(uint sz)
{_G();
    GC_ALLOCATE(NVector_O, nv );
    nv->setup(sz);
    return nv;
}




/*
 *      NVector_O
 *
 *	Author:	Christian Schafmeister (1991)
 *
 */
    NVector_O::NVector_O( const NVector_O& orig ) : NVector_O::Base(orig)
{_G();
    this->_Size = orig._Size;
    this->_Values.assign(orig._Values.begin(),orig._Values.end());
}



void	NVector_O::initialize()
{
    this->setup(0);
}

void	NVector_O::setup( uint iElements )
{_G();
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


    core::T_sp	NVector_O::deepCopy() const
{ _G();
    GC_COPY(NVector_O, vv, *this->const_sharedThis<NVector_O>() ); // = RP_Copy<NVector_O>(this->const_sharedThis<NVector_O>());
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
void	NVector_O::fill( double dValue )
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
{ _G();
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
 *      Add two vectors and return the result in (this)
 */
void	NVector_O::add( NVector_sp x, NVector_sp y )
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
void	NVector_O::XPlusYTimesScalar( NVector_sp x, NVector_sp y, double s )
{_G();
vector<double>::iterator	resi, xi,yi;
    LOG(BF("this->_Values.size() = %d") % this->_Values.size() );
    LOG(BF("x->size() = %d") % x->size() );
    LOG(BF("y->size() = %d") % y->size() );
    ASSERTP( this->_Values.size() == x->size(),"NVector_O::addTimesScalar>>mismatched size");
    ASSERT( this->_Values.size() == y->size());
    for ( resi = this->_Values.begin(), xi=x->_Values.begin(), yi=y->_Values.begin();
	    resi!=this->_Values.end(); resi++, xi++, yi++ ) {
	*resi = *xi+(*yi)*s;
    }
}

/*
 *      addScalar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a scalar to each value (this)
 */
void	NVector_O::addScalar( NVector_sp x, double s )
{_G();
vector<double>::iterator	resi, xi;
    LOG(BF("this->_Values.size() = %d") % this->_Values.size() );
    LOG(BF("x->size() = %d") % x->size() );
    ASSERTP( this->_Values.size() == x->size(),"NVector_O::addTimesScalar>>mismatched size");
    for ( resi = this->_Values.begin(), xi=x->_Values.begin();
	    resi!=this->_Values.end(); resi++, xi++ ) {
	*resi = *xi+s;
    }
}


void	NVector_O::inPlaceAddTimesScalar( NVector_sp dir, double s )
{_OF();
vector<double>::iterator	resi, diri;
    ASSERT( this->_Values.size() == dir->size());
    for ( resi = this->_Values.begin(), diri=dir->_Values.begin();
	    resi!=this->_Values.end(); resi++, diri++ ) {
	*resi += (*diri)*s;
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






/*
 *      dotProduct
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the DOT product of two vectors.
 *      return:  (this) . x
 */
double	NVector_O::dotProduct( NVector_sp x )
{
double          dDot;
vector<double>::iterator	resi,xi;
    dDot = 0.0;
    for ( resi=this->_Values.begin(),xi=x->_Values.begin(); resi!=this->_Values.end(); resi++,xi++ ) {
	dDot += (*resi)*(*xi);
    }
    return(dDot);
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
{_G();
    node->archiveVectorDouble( "vec", this->_Values );
    LOG(BF("size of this->_Values=%d") % this->_Values.size()  );
    if ( node->loading() )
    {
       this->_Size = this->_Values.size();
    }
}
#endif




/*
 *      magnitude
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the length of the vector.
 */
double	NVector_O::magnitude()
{
double  dDot;

    dDot = this->squared();
    return(sqrt(dDot));
}


double	NVector_O::rmsMagnitude()
{
double  dDot;

    dDot = this->squared();
    dDot /= (double)this->size();
    return(sqrt(dDot));
}


#if 0
/*
 * termDifferencesAsXml
 *
 * List the term-by-term differences and return as an XML object
 */
adapt::QDomNode_sp	NVector_O::termDifferencesAsXml(const string& thisName, NVector_sp n, const string& otherName)
{
    adapt::QDomNode_sp xml = adapt::QDomNode_O::create("Differences");
    adapt::QDomNode_sp part;
    uint		i;
    double	avg;
    vector<double>::iterator	resi, ni;
    resi = this->_Values.begin();
    ni = n->_Values.begin();
    for ( i=0; i<this->_Size; i++, resi++, ni++ ) {
	avg = (*resi+*ni)/2.0;
	if ( avg < 0.0000001 ) continue;
	if ( fabs(*resi-*ni)/avg > 0.1 ) {
	    part = adapt::QDomNode_O::create("V");
	    part->addAttributeInt("i",i);
	    part->addAttributeDoubleScientific(thisName,*resi);
	    part->addAttributeDoubleScientific(otherName,*ni);
	    xml->addChild(part);
	}
    }
    return xml;
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


double	NVector_O::angleWithVector(NVector_sp other)
{
#define	VERY_SMALL 1.0e-6
double	lenThis, lenOther, dot;
    lenThis = this->magnitude();
    lenOther = other->magnitude();
    dot = this->dotProduct(other);
    if ( fabs(lenThis) > VERY_SMALL && fabs(lenOther)>VERY_SMALL) {
	dot /= (lenThis*lenOther);
	if ( dot > 1.0 ) dot = 1.0;
	if ( dot < -1.0 ) dot = -1.0;
	return acos(dot);
    }
    return 0.0;
}


double	NVector_O::rmsDistanceFrom(NVector_sp v)
{_OF();
uint	i;
double	e, sum;
    ASSERT(v->size() == this->size() );
    sum = 0.0;
    for ( i=0; i<this->size(); i++ ) {
	e = this->_Values[i]-v->_Values[i];
	sum += e*e;
    }
    sum /= (float)(this->size());
    return sqrt(sum);
}


};

