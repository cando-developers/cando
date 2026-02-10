/*
    File: superposeEngine.cc
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
#define        DEBUG_LEVEL_FULL


//
//        superpose
//
//        Superpose two sets of points,
//        return the transformation that transforms moveable
//        to the fixed ones.
//




#include <stdio.h>
#include <math.h>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/geom/matrix.h>
#include <cando/geom/omatrix.h>
#include <clasp/core/symbolTable.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/geom/coordinateArray.h>
#include <cando/chem/loop.h>
#include <cando/chem/superposeEngine.h>
#include <clasp/core/wrappers.h>
#include <cando/chem/candoScript.h>
#include <cando/geom/vector3.h>

namespace chem
{

SYMBOL_EXPORT_SC_(ChemPkg,superpose);

CL_DEFUN SuperposeEngine_sp chem__make_SuperposeEngine()
{
  auto obj = gc::GC<SuperposeEngine_O>::allocate();
  return obj;
}


CL_LISPIFY_NAME("debugString");
CL_DEFMETHOD     string SuperposeEngine_O::debugString()
{
  stringstream ss;
  for (uint i=0; i<this->_FixedIndexes->length(); i++ )
  {
    int idx = (*this->_FixedIndexes)[i];
    ss << "Fixed#" << i << " " << (*this->_FixedCoordinates)[idx].asString() << std::endl;
  }
  for (uint j=0; j<this->_MoveableIndexes->length(); j++ )
  {
    int idx = (*this->_MoveableIndexes)[j];
    ss << "Moveable#" << j << " " << (*this->_MoveableCoordinates)[idx].asString() << std::endl;
  }
  return ss.str();
}

void SuperposeEngine_O::initialize()
{
    this->Base::initialize();
//    this->_FixedIndexes = core::ComplexVector_byte32_t_O::make_vector(16,0,core::clasp_make_fixnum(0));
//    this->_MoveableIndexes = core::ComplexVector_byte32_t_O::make_vector(16,0,core::clasp_make_fixnum(0));
    this->_FixedIndexes = core::ComplexVector_byte32_t_O::make_vector(16);
    this->_MoveableIndexes = core::ComplexVector_byte32_t_O::make_vector(16);
    this->_FixedCoordinates = geom::ComplexVectorCoordinate_O::make_vector(16,Vector3(),core::clasp_make_fixnum(0));
    this->_MoveableCoordinates = geom::ComplexVectorCoordinate_O::make_vector(16,Vector3(),core::clasp_make_fixnum(0));
}

void SuperposeEngine_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,FixedIndexes),this->_FixedIndexes);
  node->field(INTERN_(kw,FixedCoordinates),this->_FixedCoordinates);
  node->field(INTERN_(kw,MoveableIndexes),this->_MoveableIndexes);
  node->field(INTERN_(kw,MoveableCoordinates),this->_MoveableCoordinates);
  node->field(INTERN_(kw,Transform),this->_Transform);
}

CL_DEFMETHOD
void        SuperposeEngine_O::eraseMoveablePoints()
{
  this->_MoveableIndexes->fillPointerSet(0);
    this->_MoveableCoordinates->fillPointerSet(0);
}

CL_DEFMETHOD
void        SuperposeEngine_O::eraseFixedPoints()
{
  this->_FixedIndexes->fillPointerSet(0);
  this->_FixedCoordinates->fillPointerSet(0);
}

CL_DEFMETHOD
    void SuperposeEngine_O::appendMoveablePoint(const Vector3& pos)
{
    ASSERT(this->_MoveableIndexes.notnilp());
    ASSERT(this->_MoveableCoordinates.notnilp());
    LOG("appendMoveablePoint vector: {}" , pos.asString() );
    int idx = this->_MoveableCoordinates->length();
    this->_MoveableIndexes->vectorPushExtend(idx);
    this->_MoveableCoordinates->vectorPushExtend_Vector3(pos);
    LOG("There are now {} MoveableIndexes" , this->_MoveableIndexes->length());
    LOG("There are now {} MoveableCoordinates" , this->_MoveableCoordinates->length());
}

    void SuperposeEngine_O::appendFixedPoint(const Vector3& pos)
{
    ASSERT(this->_FixedIndexes.notnilp());
    ASSERT(this->_FixedCoordinates.notnilp());
    LOG("appendFixedPoint vector: {}" , pos.asString() );
    int idx = this->_FixedCoordinates->length();
    this->_FixedIndexes->vectorPushExtend(idx);
    this->_FixedCoordinates->vectorPushExtend_Vector3(pos);
    LOG("There are now {} FixedIndexes" , this->_FixedIndexes->length());
    LOG("There are now {} FixedCoordinates" , this->_FixedCoordinates->length());
}



CL_LISPIFY_NAME("getNumberOfFixedPoints");
CL_DEFMETHOD int        SuperposeEngine_O::getNumberOfFixedPoints()
{
    return this->_FixedIndexes->length();
}


CL_LISPIFY_NAME("getNumberOfMoveablePoints");
CL_DEFMETHOD int        SuperposeEngine_O::getNumberOfMoveablePoints()
{
    return this->_MoveableIndexes->length();
}


/*! You might look here for a reference on how the math works:
https://is.muni.cz/th/u0qtg/thesis.pdf
I found this reference in 2018 and I wrote the code in early 2010's -
I don't recall the original reference I used.
*/
CL_DEFMETHOD
void        SuperposeEngine_O::doSuperpose()
{
    VectorVector3s                        Sj;
    VectorVector3s                        Si;
    VectorVector3s::iterator        itSi,itSj;
    Vector3                                fixedCenter,vTemp;
    Vector3                                moveableCenter;
    Matrix                                mat,M, MT, trM, P,evecs,em,trans,rot;
    Vector4                                v4,evals,largestEv,quaternion;
    int                                iLargestEv;
    double                                X0, X1, X2, x0, x1, x2;
    double                                l,m,n,s,ll,mm,nn,ss;
    int                                fixedIndexesSize, moveableIndexesSize;
    fixedIndexesSize = this->_FixedIndexes->length();
    moveableIndexesSize = this->_MoveableIndexes->length();
    LOG("Moveable indexes fixed({}) moveable({})" , fixedIndexesSize , moveableIndexesSize );
    if (fixedIndexesSize != moveableIndexesSize) {
      SIMPLE_ERROR("Number of fixed points {} must match the number of moveable points {}",
                   fixedIndexesSize,
                   moveableIndexesSize );
    }
    ASSERTF(fixedIndexesSize>=3, ("You must have at least 3 points to superpose and you only have: {}") , fixedIndexesSize );
    Sj.resize(this->_FixedIndexes->length());
    if ( this->_MoveableIndexes->length() < 3 )
    {
        SIMPLE_ERROR("Number of MoveableIndexes must be greater than 3");
    }
    Si.resize(this->_MoveableIndexes->length());
    LOG("this->_FixedIndexes->length()={}" , this->_FixedIndexes->length()  );
    LOG("this->_MoveableIndexes->length()={}" , this->_MoveableIndexes->length()  );
    LOG("this->_FixedCoordinates->length() = {}\n" , this->_FixedCoordinates->length());
    {
        fixedCenter.set(0.0,0.0,0.0);
        LOG("Calculating fixedCenter\n");
        for ( size_t iaV(0); iaV<this->_FixedIndexes->length(); ++iaV) {
          LOG("  iaV = {}\n" , iaV);
          LOG("  (*this->_FixedIndexes)[iaV] = {}\n" , (*this->_FixedIndexes)[iaV] );
          LOG("  (*this->_FixedCoordinates)[(*this->_FixedIndexes)[iaV]] -> {}\n"
              , (*this->_FixedCoordinates)[(*this->_FixedIndexes)[iaV]].asString() );
          fixedCenter = (*this->_FixedCoordinates)[(*this->_FixedIndexes)[iaV]] + fixedCenter;
          LOG(" accumulating fixedCenter = {}\n" , fixedCenter.asString());
        }
        LOG("About to divide by number of indexes {}\n" , this->_FixedIndexes->length());
        fixedCenter = fixedCenter.multiplyByScalar(1.0/(double)(this->_FixedIndexes->length()));
        LOG( "Translating fixed to fixed geometric center: {} " , fixedCenter.asString() );
        VectorVector3s::iterator itS = Sj.begin();
        for ( size_t iaV(0);
              iaV<this->_FixedIndexes->length();
              iaV++,itS++ ) {
          *itS = (*this->_FixedCoordinates)[(*this->_FixedIndexes)[iaV]] - fixedCenter;
          LOG("Centered _FixedCoordinates index={}   original={}  new={}" , iaV , (*this->_FixedCoordinates)[(*this->_FixedIndexes)[iaV]].asString() , (*itS).asString() );
        }
    }
    {
        moveableCenter.set(0.0,0.0,0.0);
        for ( size_t iaV(0); iaV<this->_MoveableIndexes->length(); ++iaV) {
          moveableCenter =(*this->_MoveableCoordinates)[(*this->_MoveableIndexes)[iaV]] + moveableCenter;
        }
        moveableCenter = moveableCenter.multiplyByScalar(1.0/(double)(this->_MoveableIndexes->length()));
        LOG("Translating moveable to moveable geometric center: {}" , moveableCenter.asString().c_str()  );
        auto itS=Si.begin();
        for ( size_t iaV(0); itS!=Si.end(); iaV++,itS++ ) {
          *itS = (*this->_MoveableCoordinates)[(*this->_MoveableIndexes)[iaV]] - moveableCenter;
          LOG("Centered _MoveableCoordinates index={}   original={}  new={}" , iaV , (*this->_MoveableCoordinates)[(*this->_MoveableIndexes)[iaV]].asString() , (*itS).asString() );
        }
    }
    LOG( "fixedCenter = {}" , fixedCenter.asString() );
    LOG( "moveableCenter = {}" , moveableCenter.asString() );

    itSj = Sj.begin();
    itSi = Si.begin();
    M.setValue(0.0);
    for ( uint c = 0; c< Sj.size(); c++ ) {
        X0 = itSj->getX();
        X1 = itSj->getY();
        X2 = itSj->getZ();
        x0 = itSi->getX();
        x1 = itSi->getY();
        x2 = itSi->getZ();
        LOG("Iterating c={}" , c  );
        LOG("Iterating x0,x1,x2 ={}, {}, {}" , x0 , x1 , x2  );
        LOG("Iterating X0,X1,X2 ={}, {}, {}" , X0 , X1 , X2  );
        M.at(0, 0) += x0*X0;
        M.at(0, 1) += x0*X1;
        M.at(0, 2) += x0*X2;
        M.at(1, 0) += x1*X0;
        M.at(1, 1) += x1*X1;
        M.at(1, 2) += x1*X2;
        M.at(2, 0) += x2*X0;
        M.at(2, 1) += x2*X1;
        M.at(2, 2) += x2*X2;
        itSj++;
        itSi++;
    }
    LOG( "M= \n{}" , M.asString() );
    MT = M.transpose();
    v4.set( 0.0, 0.0, 0.0, M.trace() );
    trM.setFromQuaternion(v4);
    P = M + (MT - (trM*2.0));
    P.at(0, 3) = M.at(1, 2) - M.at(2, 1);
    P.at(3, 0) = P.at(0, 3);
    P.at(1, 3) = M.at(2, 0) - M.at(0, 2);
    P.at(3, 1) = P.at(1, 3);
    P.at(2, 3) = M.at(0, 1) -M.at(1, 0);
    P.at(3, 2) = P.at(2, 3);
    P.at(3, 3) = 0.0;
    LOG("P=\n{}" , P.asString() );
    P.eigenSystem( evals, evecs );
    iLargestEv = evals.indexOfLargestElement();
    largestEv = evecs.getCol(iLargestEv);
    LOG( "largestEv=\n{}" , largestEv.asString() );
    em.setFromQuaternion(largestEv);

    LOG( "eigenMatrix=\n{}" , em.asString() );
    quaternion = em.getCol(3);
    LOG( "quaternion = \n{}" , quaternion.asString() );
    //
    // Convert the quaternion into a rotation matrix
    rot.setValue(0.0);
    l = quaternion.getW();
    m = quaternion.getX();
    n = quaternion.getY();
    s = quaternion.getZ();
    ll = l*l;
    mm = m*m;
    nn = n*n;
    ss = s*s;
    rot.at(0, 0) = ll - mm - nn + ss;
    rot.at(0, 1) = 2*(l*m - n*s);
    rot.at(0, 2) = 2*(l*n + m*s);
    rot.at(1, 0) = 2*(l*m + n*s);
    rot.at(1, 1) = -ll + mm - nn + ss;
    rot.at(1, 2) = 2*(m*n - l*s);
    rot.at(2, 0) = 2*(l*n - m*s);
    rot.at(2, 1) = 2*(m*n + l*s);
    rot.at(2, 2) = -ll - mm + nn + ss;
    rot.at(3, 3) = 1.0;
    LOG( "rot=\n{}" , rot.asString());
    vTemp = moveableCenter*-1.0;
    trans.translate(vTemp);
    LOG( "moveableTrans=\n{}" , trans.asString() );
//    mat.setFromQuaternion(quaternion);
    mat = rot*trans;
    LOG( "mat*moveableTrans= \n{}" , mat.asString() );
    trans.translate(fixedCenter);
    LOG( "fixedTrans = \n{}" , trans.asString() );
    this->_Transform = trans*mat;
    LOG("this->_Transform = \n{}" , this->_Transform.asString());
}


#if 0
/*!
  Return the rootMeanSquare difference between the two
  collections of points.
*/
double        SuperposeEngine_O::sumOfSquaresOfDifferences(ScorerState_sp scorerState )
{
    core::ComplexVector_byte32_t_O::iterator                itFixed;
    core::ComplexVector_byte32_t_O::iterator                ititMoved;
    Vector3                                moved, diff;
    double                                sum;

    ASSERTNOTNULL(this->_MoveableIndexes);
    ASSERTNOTNULL(this->_FixedIndexes);
    LOG_SCORE(scorerState,("SuperposeEngine_O::sumOfSquaresOfDifferences{"));
#if        DEBUG_SCORE_EVALUATION
    LOG_SCORE(scorerState,("There are %d moveable indexes") , this->_MoveableIndexes->length() );
    LOG_SCORE(scorerState,("There are %d fixed indexes") , this->_FixedIndexes->length() );
    for ( ititMoved = this->_MoveableIndexes->begin(),
              itFixed = this->_FixedIndexes->begin();
          itFixed != this->_FixedIndexes->end();
          ititMoved++, itFixed++ )
    {
      LOG_SCORE(scorerState,"Superpose fixed%s with moveable%s" ,
                this->_FixedCoordinates->getElement(*itFixed).asString()
                , this->_MoveableCoordinates->getElement(*ititMoved).asString()
            );
    }
    LOG_SCORE(BF(scorerState,"SuperposeEngine_O calculated transform\n%s")% this->_Transform.asString());
#endif
    double result = this->sumOfSquaresOfDifferences();
    LOG_SCORE(scorerState,BF("SuperposeEngine_O::returning result=%lf")% result);
    LOG_SCORE(scorerState,BF("SuperposeEngine_O::sumOfSquaresOfDifferences}"));
    return result;
}
#endif

/*!
  Return the rootMeanSquare difference between the two
  collections of points.
*/
CL_DEFMETHOD
double        SuperposeEngine_O::sumOfSquaresOfDifferences()
{
  Vector3                                moved, diff;
  double                                sum;

    LOG("Dump of coordinates being that I will superpose" );
    ASSERTNOTNULL(this->_MoveableIndexes);
    ASSERTNOTNULL(this->_FixedIndexes);
#ifdef        DEBUG_ON
    LOG("There are {} moveable indexes" , this->_MoveableIndexes->length()  );
    LOG("There are {} fixed indexes" , this->_FixedIndexes->length()  );
    size_t ititMoved;
    size_t itFixed;
    for ( ititMoved=0, itFixed=0;
          itFixed < this->_FixedIndexes->length();
          ititMoved++, itFixed++ )
    {
        LOG("Superpose fixed{} with moveable{}"
            , (*this->_FixedCoordinates)[itFixed].asString()
            , (*this->_MoveableCoordinates)[ititMoved].asString() );
    }
#endif
    LOG("Superpose transform: \n{}\n" , this->_Transform.asString());
    ASSERT(this->_FixedIndexes->length()!=0);
    sum = 0.0;
    for ( size_t ititMoved(0), itFixed(0);
          itFixed < this->_FixedIndexes->length();
          ititMoved++, itFixed++ ) {
      moved = this->_Transform.multiplyByVector3((*this->_MoveableCoordinates)[ititMoved]);
      diff = (*this->_FixedCoordinates)[itFixed]-moved;
      LOG("Looking at itFixed {} {}   ititMoved {} {}  -->  diff -> {}\n" , itFixed , (*this->_FixedCoordinates)[itFixed].asString() , ititMoved , moved.asString() , diff.asString());
      sum = sum + diff.dotProduct(diff);
    }
    LOG("Calculated sum={}" , sum  );
    return sum;
}

//
//        rootMeanSquareDifference
//
//        Return the rootMeanSquare difference between the two
//        collections of points.
CL_LISPIFY_NAME("rootMeanSquareDifference");
CL_DEFMETHOD double        SuperposeEngine_O::rootMeanSquareDifference()
{
    double sumOfSquares = this->sumOfSquaresOfDifferences();
    LOG("Number of moveable indexes = {}" , this->_MoveableIndexes->length() );
    ASSERT_gt(this->_MoveableIndexes->length(),0);
    LOG("sumOfSquares = {}" , sumOfSquares );
    double meanSumOfSquares = sumOfSquares/this->_MoveableIndexes->length();
    LOG("meanSumOfSquares = {}" , meanSumOfSquares );
    double rms = sqrt(meanSumOfSquares);
    LOG("rms = {}" , rms );
    return rms;
}

CL_LISPIFY_NAME("setFixedPoints");
CL_DEFMETHOD void SuperposeEngine_O::setFixedPoints(core::ComplexVector_byte32_t_sp fi, geom::SimpleVectorCoordinate_sp fc) {
    LOG("SuperposeEngine_O::setFixedPoints --> number of points={}", fc->length());
    if (!((fi->length() >= 3) && ((fi->length() <= fc->length())))) {
      SIMPLE_ERROR("In SuperposeEngine setFixedPoints - number of indexes fi->length() = {} must be at least 3\nNumber of "
                   "coordinates --> fc->length() = {} must be >= to fi->size({})",
                   fi->length(), fc->length(), fi->length());
    }
    this->_FixedIndexes = gc::As<core::ComplexVector_byte32_t_sp>(fi);
    geom::SimpleVectorCoordinate_sp data = fc->copy();
    this->_FixedCoordinates = geom::ComplexVectorCoordinate_O::make_vector(data);
}

CL_LISPIFY_NAME("setFixedAllPoints");
CL_DEFMETHOD void        SuperposeEngine_O::setFixedAllPoints( core::Array_sp afc )
{
  if (gc::IsA<geom::SimpleVectorCoordinate_sp>(afc) ) {
    auto fc = gc::As_unsafe<geom::SimpleVectorCoordinate_sp>(afc);
    this->_FixedCoordinates = geom::ComplexVectorCoordinate_O::make_vector(fc->length(),Vector3(),core::make_fixnum(fc->length()));
    this->_FixedIndexes = core::ComplexVector_byte32_t_O::make_vector(fc->length());
    size_t ii;
    for ( ii=0; ii<fc->length(); ii++ ) {
      (*this->_FixedCoordinates)[ii] = (*fc)[ii];
      (*this->_FixedIndexes)[ii] = ii;
//    printf("%s:%d Set fixed coordinate[%zu] -> %lf, %lf, %lf\n", __FILE__, __LINE__, ii, (*fc)[ii].getX(), (*fc)[ii].getY(), (*fc)[ii].getZ());
      LOG("setFixedAllPoints index@{} " , ii  );
    }
    return;
  } else if (gc::IsA<NVector_sp>(afc) ) {
    auto fc = gc::As_unsafe<NVector_sp>(afc);
    ASSERTF(fc->length()>=3*3,("You must have at least three moveable points and there are only {}") , fc->length()/3 );
    this->_FixedCoordinates = geom::ComplexVectorCoordinate_O::make_vector(fc->length()/3,Vector3(),core::make_fixnum(fc->length()/3));
    this->_FixedIndexes = core::ComplexVector_byte32_t_O::make_vector(fc->length()/3);
    uint                ii;
    for ( ii=0; ii<fc->length()/3; ii++ )
    {
      (*this->_FixedCoordinates)[ii].set((*fc)[(ii*3)+0], (*fc)[(ii*3)+1], (*fc)[(ii*3)+2] );
      (*this->_FixedIndexes)[ii] = ii;
      LOG("setMoveableAllPoints index@{} " , ii  );
//    printf("%s:%d Set moveable coordinate[%d] -> %lf, %lf, %lf\n", __FILE__, __LINE__, ii, (*fc)[ii].getX(), (*fc)[ii].getY(), (*fc)[ii].getZ());
    }
    return;
  }
  TYPE_ERROR(afc,core::Cons_O::createList(NVector_O::static_classSymbol(), geom::SimpleVectorCoordinate_O::static_classSymbol()));
}



CL_LISPIFY_NAME("setMoveablePoints");
CL_DEFMETHOD void        SuperposeEngine_O::setMoveablePoints( core::ComplexVector_byte32_t_sp mi, geom::SimpleVectorCoordinate_sp mc )
{
    LOG("SuperposeEngine_O::setMoveablePoints --> number of points={}" , mc->length()  );
    ASSERTP( (mi->length()>=3), "There must be at least three indexes" );
    ASSERTP( mi->length()<=mc->length(), "There must be at least as many coordinates as indexes");
    this->_MoveableIndexes = gc::As<core::ComplexVector_byte32_t_sp>(mi);
    geom::SimpleVectorCoordinate_sp data = mc->copy();
    this->_MoveableCoordinates = geom::ComplexVectorCoordinate_O::make_vector(data);
}

CL_LISPIFY_NAME("setMoveableAllPoints");
CL_DEFMETHOD void        SuperposeEngine_O::setMoveableAllPoints( core::Array_sp amc )
{
  if (gc::IsA<geom::SimpleVectorCoordinate_sp>(amc) ) {
    auto mc = gc::As_unsafe<geom::SimpleVectorCoordinate_sp>(amc);
    uint                ii;
    ASSERTF(mc->length()>=3,("You must have at least three moveable points and there are only {}") , mc->length() );
    this->_MoveableCoordinates = geom::ComplexVectorCoordinate_O::make_vector(mc->length(),Vector3(),core::make_fixnum(mc->length()));
    this->_MoveableIndexes = core::ComplexVector_byte32_t_O::make_vector(mc->length());
    for ( ii=0; ii<mc->length(); ii++ )
    {
      (*this->_MoveableCoordinates)[ii] = (*mc)[ii];
      (*this->_MoveableIndexes)[ii] = ii;
      LOG("setMoveableAllPoints index@{} " , ii  );
//    printf("%s:%d Set moveable coordinate[%d] -> %lf, %lf, %lf\n", __FILE__, __LINE__, ii, (*mc)[ii].getX(), (*mc)[ii].getY(), (*mc)[ii].getZ());
    }
    return;
  } else if (gc::IsA<NVector_sp>(amc)) {
    auto mc = gc::As_unsafe<NVector_sp>(amc);
    ASSERTF(mc->length()>=3*3,("You must have at least three moveable points and there are only {}") , mc->length() );
    this->_MoveableCoordinates = geom::ComplexVectorCoordinate_O::make_vector(mc->length()/3,Vector3(),core::make_fixnum(mc->length()));
    this->_MoveableIndexes = core::ComplexVector_byte32_t_O::make_vector(mc->length()/3);
    uint                ii;
    for ( ii=0; ii<mc->length()/3; ii++ )
    {
      (*this->_MoveableCoordinates)[ii].set((*mc)[(ii*3)+0], (*mc)[(ii*3)+1], (*mc)[(ii*3)+2] );
      (*this->_MoveableIndexes)[ii] = ii;
      LOG("setMoveableAllPoints index@{} " , ii  );
//    printf("%s:%d Set moveable coordinate[%d] -> %lf, %lf, %lf\n", __FILE__, __LINE__, ii, (*mc)[ii].getX(), (*mc)[ii].getY(), (*mc)[ii].getZ());
    }
    return;
  }
  TYPE_ERROR(amc,core::Cons_O::createList(NVector_O::static_classSymbol(), geom::SimpleVectorCoordinate_O::static_classSymbol()));
}



CL_LISPIFY_NAME("superpose");
CL_DEFMETHOD Matrix        SuperposeEngine_O::superpose()
{
    LOG("SuperposeEngine_O::superpose()" );
    this->doSuperpose();
    LOG("Carried out superpose and the transform is:{}" , (this->_Transform.asString().c_str() ) );
    return this->_Transform;
}








void SuperposeSelectedAtoms_O::initialize()
{
    this->Base::initialize();
    this->_Matter = nil<Matter_O>();
    this->_SuperposeAtoms.clear();
}


#ifdef XML_ARCHIVE
void SuperposeSelectedAtoms_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("matter",this->_Matter);
    node->archiveVector0("superposeAtoms",this->_SuperposeAtoms);
}
#endif


void SuperposeSelectedAtoms_O::updateSuperposeAtoms()
{
    Loop lAtoms;
    this->_SuperposeAtoms.clear();
    lAtoms.loopTopGoal(this->_Matter,ATOMS);
    uint numAtoms = 0;
    uint numSelectedAtoms = 0;
    SYMBOL_EXPORT_SC_(ChemPkg,superpose);
    core::Symbol_sp superposeSymbol = _sym_superpose;
    while ( lAtoms.advance() )
    {
        Atom_sp a = lAtoms.getAtom();
        if ( a->getPropertyOrDefault(superposeSymbol,_lisp->_true()).isTrue() )
        {
            this->_SuperposeAtoms.push_back(a);
            numSelectedAtoms++;
        }
        numAtoms++;
    }
    ASSERT_gt(numAtoms,0);
    ASSERT_gt(numSelectedAtoms,0);
    ASSERT_gt(this->_SuperposeAtoms.size(),0);
}

CL_LISPIFY_NAME("setSuperposeMatter");
CL_DEFMETHOD void SuperposeSelectedAtoms_O::setMatter(Matter_sp matter)
{
  this->_Matter = matter->copy(nil<core::T_O>());
  this->_SuperposeAtoms.clear();
}

CL_LISPIFY_NAME("setMatterWithSelectedAtoms");
CL_DEFMETHOD void SuperposeSelectedAtoms_O::setMatterWithSelectedAtoms(Matter_sp matter)
{
    this->setMatter(matter);
    this->updateSuperposeAtoms();
}


CL_LISPIFY_NAME("extractCoordinates");
CL_DEFMETHOD geom::SimpleVectorCoordinate_sp SuperposeSelectedAtoms_O::extractCoordinates(Matter_sp matter)
{
    if ( !this->_Matter->equal(matter) )
    {
        SIMPLE_ERROR("The Matters are not equal");
    }
    this->_Matter->transferCoordinates(matter);
    ASSERT_gt(this->_SuperposeAtoms.size(),0);
    geom::SimpleVectorCoordinate_sp coords = geom::SimpleVectorCoordinate_O::make(this->_SuperposeAtoms.size());
    uint idx = 0;
    for (gctools::Vec0<Atom_sp>::iterator ai=this->_SuperposeAtoms.begin(); ai!=this->_SuperposeAtoms.end(); ai++, idx++ )
    {
      (*coords)[idx] = (*ai)->getPosition();
    }
    return coords;
}

CL_LISPIFY_NAME("getSuperposeMatter");
CL_DEFMETHOD Matter_sp SuperposeSelectedAtoms_O::getMatter()
{
    return this->_Matter;
}

CL_LISPIFY_NAME("copyMatterCoordinatesIntoFixedCoordinates");
CL_DEFMETHOD void SuperposeSelectedAtoms_O::copyMatterCoordinatesIntoFixedCoordinates(Matter_sp matter)
{
    geom::SimpleVectorCoordinate_sp ca = this->extractCoordinates(matter);
    ASSERT_gt(ca->length(),0);
    this->setFixedAllPoints(ca);
}

CL_LISPIFY_NAME("copyMatterCoordinatesIntoMoveableCoordinates");
CL_DEFMETHOD void SuperposeSelectedAtoms_O::copyMatterCoordinatesIntoMoveableCoordinates(Matter_sp matter)
{
    geom::SimpleVectorCoordinate_sp ca = this->extractCoordinates(matter);
    ASSERT_gt(ca->length(),0);
    this->setMoveableAllPoints(ca);
}


DOCGROUP(cando);
CL_DEFUN void chem__align_matter(Matter_sp moveable_matter, core::T_sp moveable_atoms, core::T_sp fixed_atoms)
{
  size_t num_moveable = core::cl__length(moveable_atoms);
  size_t num_fixed = core::cl__length(fixed_atoms);
  if (num_moveable<3) {
    SIMPLE_ERROR("There must be at least 3 moveable atoms - {} were passed" , num_moveable);
  }
  if (num_fixed<3) {
    SIMPLE_ERROR("There must be at least 3 fixed atoms - {} were passed" , num_fixed);
  }
  if (num_moveable != num_fixed) {
    SIMPLE_ERROR("The number of moveable points {} must match the number of fixed points" , num_moveable , num_fixed);
  }
  geom::SimpleVectorCoordinate_sp coords_fixed = chem__make_simple_vector_coordinate_from_atom_sequence(fixed_atoms);
  geom::SimpleVectorCoordinate_sp coords_moveable = chem__make_simple_vector_coordinate_from_atom_sequence(moveable_atoms);
  SuperposeEngine_sp engine = SuperposeEngine_O::create();
  engine->setFixedAllPoints(coords_fixed);
  engine->setMoveableAllPoints(coords_moveable);
  Matrix transform = engine->superpose();
  moveable_matter->applyTransformToAtoms(transform);
}

DOCGROUP(cando);
CL_DEFUN Matrix chem__align_transform(core::T_sp moveable_atoms, core::T_sp fixed_atoms)
{
  size_t num_moveable = core::cl__length(moveable_atoms);
  size_t num_fixed = core::cl__length(fixed_atoms);
  if (num_moveable<3) {
    SIMPLE_ERROR("There must be at least 3 moveable atoms - {} were passed" , num_moveable);
  }
  if (num_fixed<3) {
    SIMPLE_ERROR("There must be at least 3 fixed atoms - {} were passed" , num_fixed);
  }
  if (num_moveable != num_fixed) {
    SIMPLE_ERROR("The number of moveable points {} must match the number of fixed points" , num_moveable , num_fixed);
  }
  geom::SimpleVectorCoordinate_sp coords_fixed = chem__make_simple_vector_coordinate_from_atom_sequence(fixed_atoms);
  geom::SimpleVectorCoordinate_sp coords_moveable = chem__make_simple_vector_coordinate_from_atom_sequence(moveable_atoms);
  SuperposeEngine_sp engine = SuperposeEngine_O::create();
  engine->setFixedAllPoints(coords_fixed);
  engine->setMoveableAllPoints(coords_moveable);
  Matrix transform = engine->superpose();
  return transform;
}

#if 0
#define SUPERPOSE_LOG(x) core::lisp_write(x)
#else
#define SUPERPOSE_LOG(x)
#endif

CL_DOCSTRING("Return a matrix that best aligns the moveable-coordinates with the fixed-coordinates using the indexes3 atoms");
CL_DEFUN
Matrix chem__align_coordinates(core::SimpleVector_byte32_t_sp indexes3, NVector_sp fixedCoordinates, NVector_sp moveableCoordinates )
{
  VectorVector3s                        Sj;
  VectorVector3s                        Si;
  VectorVector3s::iterator        itSi,itSj;
  Vector3                                fixedCenter,vTemp;
  Vector3                                moveableCenter;
  Matrix                                mat,M, MT, trM, P,evecs,em,trans,rot;
  Vector4                                v4,evals,largestEv,quaternion;
  int                                iLargestEv;
  double                                X0, X1, X2, x0, x1, x2;
  double                                l,m,n,s,ll,mm,nn,ss;
  int                                fixedIndexesSize, moveableIndexesSize;

  Sj.resize(indexes3->length());
  if ( indexes3->length() < 3 ) {
    SIMPLE_ERROR("Number of indexes3 must be greater than 3");
  }
  Si.resize(indexes3->length());
  {
    fixedCenter.set(0.0,0.0,0.0);
    for ( size_t iaV(0); iaV<indexes3->length(); ++iaV) {
      Vector3 tv;
      tv.set((*fixedCoordinates)[(*indexes3)[iaV]+0],
             (*fixedCoordinates)[(*indexes3)[iaV]+1],
             (*fixedCoordinates)[(*indexes3)[iaV]+2]);
      fixedCenter = fixedCenter + tv;
    }
    fixedCenter = fixedCenter.multiplyByScalar(1.0/(double)(indexes3->length()));
    VectorVector3s::iterator itS = Sj.begin();
    for ( size_t iaV(0); iaV<indexes3->length(); iaV++,itS++ ) {
      Vector3 tv;
      tv.set((*fixedCoordinates)[(*indexes3)[iaV]+0],
             (*fixedCoordinates)[(*indexes3)[iaV]+1],
             (*fixedCoordinates)[(*indexes3)[iaV]+2]);
      *itS = tv - fixedCenter;
    }
  }
  {
    moveableCenter.set(0.0,0.0,0.0);
    for ( size_t iaV(0); iaV<indexes3->length(); ++iaV) {
      Vector3 tv;
      tv.set((*moveableCoordinates)[(*indexes3)[iaV]+0],
             (*moveableCoordinates)[(*indexes3)[iaV]+1],
             (*moveableCoordinates)[(*indexes3)[iaV]+2]);
      moveableCenter = moveableCenter + tv;
    }
    moveableCenter = moveableCenter.multiplyByScalar(1.0/(double)(indexes3->length()));
    VectorVector3s::iterator itS = Si.begin();
    for ( size_t iaV(0); iaV<indexes3->length(); iaV++,itS++ ) {
      Vector3 tv;
      tv.set((*moveableCoordinates)[(*indexes3)[iaV]+0],
             (*moveableCoordinates)[(*indexes3)[iaV]+1],
             (*moveableCoordinates)[(*indexes3)[iaV]+2]);
      SUPERPOSE_LOG(fmt::format("fixedCoordinates[{}]    = {}\n" , iaV, tv.asString() ));
      SUPERPOSE_LOG(fmt::format("moveableCoordinates[{}] = {}\n" , iaV, tv.asString() ));
      *itS = tv - moveableCenter;
    }
  }
  itSi = Si.begin();
  itSj = Sj.begin();
  M.setValue(0.0);
  for ( uint c = 0; c< Sj.size(); c++ ) {
    X0 = itSj->getX();
    X1 = itSj->getY();
    X2 = itSj->getZ();
    x0 = itSi->getX();
    x1 = itSi->getY();
    x2 = itSi->getZ();
    SUPERPOSE_LOG(fmt::format("Iterating c={}\n" , c  ));
    SUPERPOSE_LOG(fmt::format("Iterating x0,x1,x2 ={}, {}, {}\n" , x0 , x1 , x2  ));
    SUPERPOSE_LOG(fmt::format("Iterating X0,X1,X2 ={}, {}, {}\n" , X0 , X1 , X2  ));
    M.at(0, 0) += x0*X0;
    M.at(0, 1) += x0*X1;
    M.at(0, 2) += x0*X2;
    M.at(1, 0) += x1*X0;
    M.at(1, 1) += x1*X1;
    M.at(1, 2) += x1*X2;
    M.at(2, 0) += x2*X0;
    M.at(2, 1) += x2*X1;
    M.at(2, 2) += x2*X2;
    itSj++;
    itSi++;
  }
  SUPERPOSE_LOG(fmt::format( "M= {}\n" , M.asString() ));
  MT = M.transpose();
  v4.set( 0.0, 0.0, 0.0, M.trace() );
  trM.setFromQuaternion(v4);
  P = M + (MT - (trM*2.0));
  P.at(0, 3) = M.at(1, 2) - M.at(2, 1);
  P.at(3, 0) = P.at(0, 3);
  P.at(1, 3) = M.at(2, 0) - M.at(0, 2);
  P.at(3, 1) = P.at(1, 3);
  P.at(2, 3) = M.at(0, 1) -M.at(1, 0);
  P.at(3, 2) = P.at(2, 3);
  P.at(3, 3) = 0.0;
  SUPERPOSE_LOG(fmt::format("P={}\n" , P.asString() ));

  P.eigenSystem( evals, evecs );
  iLargestEv = evals.indexOfLargestElement();
  largestEv = evecs.getCol(iLargestEv);
  SUPERPOSE_LOG(fmt::format("largestEv={}\n" , largestEv.asString() ));
  em.setFromQuaternion(largestEv);

  SUPERPOSE_LOG(fmt::format( "eigenMatrix={}\n" , em.asString() ));
  quaternion = em.getCol(3);
  SUPERPOSE_LOG(fmt::format( "quaternion = {}\n" , quaternion.asString() ));
  //
  // Convert the quaternion into a rotation matrix
  rot.setValue(0.0);
  l = quaternion.getW();
  m = quaternion.getX();
  n = quaternion.getY();
  s = quaternion.getZ();
  SUPERPOSE_LOG(fmt::format("quat= {:.4f}, {:.4f}, {:.4f}, {:.4f} \n" , l, m, n, s ));
  ll = l*l;
  mm = m*m;
  nn = n*n;
  ss = s*s;
  rot.at(0, 0) = ll - mm - nn + ss;
  rot.at(0, 1) = 2*(l*m - n*s);
  rot.at(0, 2) = 2*(l*n + m*s);
  rot.at(1, 0) = 2*(l*m + n*s);
  rot.at(1, 1) = -ll + mm - nn + ss;
  rot.at(1, 2) = 2*(m*n - l*s);
  rot.at(2, 0) = 2*(l*n - m*s);
  rot.at(2, 1) = 2*(m*n + l*s);
  rot.at(2, 2) = -ll - mm + nn + ss;
  rot.at(3, 3) = 1.0;
  SUPERPOSE_LOG(fmt::format("rot=\n{}\n" , rot.asString()));
  vTemp = moveableCenter*-1.0;
  trans.translate(vTemp);
  SUPERPOSE_LOG(fmt::format( "moveableTrans={}\n" , trans.asString() ));
  //    mat.setFromQuaternion(quaternion);
  mat = rot*trans;
  SUPERPOSE_LOG(fmt::format( "mat*moveableTrans= {}\n" , mat.asString() ));
  trans.translate(fixedCenter);
  SUPERPOSE_LOG(fmt::format( "fixedTrans = {}\n" , trans.asString() ));
  Matrix transform = trans*mat;
  return transform;
}

CL_DEFUN double chem__root_mean_square_difference_transformed(core::SimpleVector_byte32_t_sp indexes3, NVector_sp fixedCoordinates, NVector_sp moveableCoordinates, Matrix transform )
{
  double sum = 0.0;
  for ( size_t iaV(0); iaV<indexes3->length(); ++iaV) {
    Vector3 fv;
    fv.set((*fixedCoordinates)[(*indexes3)[iaV]+0],
           (*fixedCoordinates)[(*indexes3)[iaV]+1],
           (*fixedCoordinates)[(*indexes3)[iaV]+2]);
    Vector3 mv;
    mv.set((*moveableCoordinates)[(*indexes3)[iaV]+0],
           (*moveableCoordinates)[(*indexes3)[iaV]+1],
           (*moveableCoordinates)[(*indexes3)[iaV]+2]);
    mv = transform.multiplyByVector3(mv);
    Vector3 diff = fv-mv;
    sum = sum + diff.dotProduct(diff);
  }
  double meanSumOfSquares = sum/indexes3->length();
  return sqrt(meanSumOfSquares);
}

};
