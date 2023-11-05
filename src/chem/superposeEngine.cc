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
  for (uint i=0; i<this->_FixedIndices->length(); i++ )
  {
    int idx = (*this->_FixedIndices)[i];
    ss << "Fixed#" << i << " " << (*this->_FixedCoordinates)[idx].asString() << std::endl;
  }
  for (uint j=0; j<this->_MoveableIndices->length(); j++ )
  {
    int idx = (*this->_MoveableIndices)[j];
    ss << "Moveable#" << j << " " << (*this->_MoveableCoordinates)[idx].asString() << std::endl;
  }
  return ss.str();
}

void SuperposeEngine_O::initialize()
{
    this->Base::initialize();
//    this->_FixedIndices = core::ComplexVector_byte32_t_O::make_vector(16,0,core::clasp_make_fixnum(0));
//    this->_MoveableIndices = core::ComplexVector_byte32_t_O::make_vector(16,0,core::clasp_make_fixnum(0));
    this->_FixedIndices = core::ComplexVector_byte32_t_O::make_vector(16);
    this->_MoveableIndices = core::ComplexVector_byte32_t_O::make_vector(16);
    this->_FixedCoordinates = geom::ComplexVectorCoordinate_O::make_vector(16,Vector3(),core::clasp_make_fixnum(0));
    this->_MoveableCoordinates = geom::ComplexVectorCoordinate_O::make_vector(16,Vector3(),core::clasp_make_fixnum(0));
}

void SuperposeEngine_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,FixedIndices),this->_FixedIndices);
  node->field(INTERN_(kw,FixedCoordinates),this->_FixedCoordinates);
  node->field(INTERN_(kw,MoveableIndices),this->_MoveableIndices);
  node->field(INTERN_(kw,MoveableCoordinates),this->_MoveableCoordinates);
  node->field(INTERN_(kw,Transform),this->_Transform);
}

CL_DEFMETHOD
void        SuperposeEngine_O::eraseMoveablePoints()
{
  this->_MoveableIndices->fillPointerSet(0);
    this->_MoveableCoordinates->fillPointerSet(0);
}

CL_DEFMETHOD
void        SuperposeEngine_O::eraseFixedPoints()
{
  this->_FixedIndices->fillPointerSet(0);
  this->_FixedCoordinates->fillPointerSet(0);
}

CL_DEFMETHOD
    void SuperposeEngine_O::appendMoveablePoint(const Vector3& pos)
{
    ASSERT(this->_MoveableIndices.notnilp());
    ASSERT(this->_MoveableCoordinates.notnilp());
    LOG("appendMoveablePoint vector: {}" , pos.asString() );
    int idx = this->_MoveableCoordinates->length();
    this->_MoveableIndices->vectorPushExtend(idx);
    this->_MoveableCoordinates->vectorPushExtend_Vector3(pos);
    LOG("There are now {} MoveableIndices" , this->_MoveableIndices->length());
    LOG("There are now {} MoveableCoordinates" , this->_MoveableCoordinates->length());
}

    void SuperposeEngine_O::appendFixedPoint(const Vector3& pos)
{
    ASSERT(this->_FixedIndices.notnilp());
    ASSERT(this->_FixedCoordinates.notnilp());
    LOG("appendFixedPoint vector: {}" , pos.asString() );
    int idx = this->_FixedCoordinates->length();
    this->_FixedIndices->vectorPushExtend(idx);
    this->_FixedCoordinates->vectorPushExtend_Vector3(pos);
    LOG("There are now {} FixedIndices" , this->_FixedIndices->length());
    LOG("There are now {} FixedCoordinates" , this->_FixedCoordinates->length());
}



CL_LISPIFY_NAME("getNumberOfFixedPoints");
CL_DEFMETHOD int        SuperposeEngine_O::getNumberOfFixedPoints()
{
    return this->_FixedIndices->length();
}


CL_LISPIFY_NAME("getNumberOfMoveablePoints");
CL_DEFMETHOD int        SuperposeEngine_O::getNumberOfMoveablePoints()
{
    return this->_MoveableIndices->length();
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
    int                                fixedIndicesSize, moveableIndicesSize;
    fixedIndicesSize = this->_FixedIndices->length();
    moveableIndicesSize = this->_MoveableIndices->length();
    LOG("Moveable indices fixed({}) moveable({})" , fixedIndicesSize , moveableIndicesSize );
    if (fixedIndicesSize != moveableIndicesSize) {
      SIMPLE_ERROR("Number of fixed points {} must match the number of moveable points {}",
                   fixedIndicesSize,
                   moveableIndicesSize );
    }
    ASSERTF(fixedIndicesSize>=3, ("You must have at least 3 points to superpose and you only have: {}") , fixedIndicesSize );
    Sj.resize(this->_FixedIndices->length());
    if ( this->_MoveableIndices->length() < 3 )
    {
        SIMPLE_ERROR("Number of MoveableIndices must be greater than 3");
    }
    Si.resize(this->_MoveableIndices->length());
    LOG("this->_FixedIndices->length()={}" , this->_FixedIndices->length()  );
    LOG("this->_MoveableIndices->length()={}" , this->_MoveableIndices->length()  );
    LOG("this->_FixedCoordinates->length() = {}\n" , this->_FixedCoordinates->length());
    {
        fixedCenter.set(0.0,0.0,0.0);
        LOG("Calculating fixedCenter\n");
        for ( size_t iaV(0); iaV<this->_FixedIndices->length(); ++iaV) {
          LOG("  iaV = {}\n" , iaV);
          LOG("  (*this->_FixedIndices)[iaV] = {}\n" , (*this->_FixedIndices)[iaV] );
          LOG("  (*this->_FixedCoordinates)[(*this->_FixedIndices)[iaV]] -> {}\n"
              , (*this->_FixedCoordinates)[(*this->_FixedIndices)[iaV]].asString() );
          fixedCenter = (*this->_FixedCoordinates)[(*this->_FixedIndices)[iaV]] + fixedCenter;
          LOG(" accumulating fixedCenter = {}\n" , fixedCenter.asString());
        }
        LOG("About to divide by number of indices {}\n" , this->_FixedIndices->length());
        fixedCenter = fixedCenter.multiplyByScalar(1.0/(double)(this->_FixedIndices->length()));
        LOG( "Translating fixed to fixed geometric center: {} " , fixedCenter.asString() );
        VectorVector3s::iterator itS = Sj.begin();
        for ( size_t iaV(0);
              iaV<this->_FixedIndices->length();
              iaV++,itS++ ) {
          *itS = (*this->_FixedCoordinates)[(*this->_FixedIndices)[iaV]] - fixedCenter;
          LOG("Centered _FixedCoordinates index={}   original={}  new={}" , iaV , (*this->_FixedCoordinates)[(*this->_FixedIndices)[iaV]].asString() , (*itS).asString() );
        }
    }
    {
        moveableCenter.set(0.0,0.0,0.0);
        for ( size_t iaV(0); iaV<this->_MoveableIndices->length(); ++iaV) {
          moveableCenter =(*this->_MoveableCoordinates)[(*this->_MoveableIndices)[iaV]] + moveableCenter;
        }
        moveableCenter = moveableCenter.multiplyByScalar(1.0/(double)(this->_MoveableIndices->length()));
        LOG("Translating moveable to moveable geometric center: {}" , moveableCenter.asString().c_str()  );
        auto itS=Si.begin();
        for ( size_t iaV(0); itS!=Si.end(); iaV++,itS++ ) {
          *itS = (*this->_MoveableCoordinates)[(*this->_MoveableIndices)[iaV]] - moveableCenter;
          LOG("Centered _MoveableCoordinates index={}   original={}  new={}" , iaV , (*this->_MoveableCoordinates)[(*this->_MoveableIndices)[iaV]].asString() , (*itS).asString() );
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

    ASSERTNOTNULL(this->_MoveableIndices);
    ASSERTNOTNULL(this->_FixedIndices);
    LOG_SCORE(scorerState,("SuperposeEngine_O::sumOfSquaresOfDifferences{"));
#if        DEBUG_SCORE_EVALUATION
    LOG_SCORE(scorerState,("There are %d moveable indices") , this->_MoveableIndices->length() );
    LOG_SCORE(scorerState,("There are %d fixed indices") , this->_FixedIndices->length() );
    for ( ititMoved = this->_MoveableIndices->begin(),
              itFixed = this->_FixedIndices->begin();
          itFixed != this->_FixedIndices->end();
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
    ASSERTNOTNULL(this->_MoveableIndices);
    ASSERTNOTNULL(this->_FixedIndices);
#ifdef        DEBUG_ON
    LOG("There are {} moveable indices" , this->_MoveableIndices->length()  );
    LOG("There are {} fixed indices" , this->_FixedIndices->length()  );
    for ( ititMoved=0, itFixed=0;
          itFixed < this->_FixedIndices->length();
          ititMoved++, itFixed++ )
    {
        LOG("Superpose fixed{} with moveable{}"
            , (*this->_FixedCoordinates)[itFixed].asString()
            , (*this->_MoveableCoordinates)[ititMoved].asString() );
    }
#endif
    LOG("Superpose transform: \n{}\n" , this->_Transform.asString());
    ASSERT(this->_FixedIndices->length()!=0);
    sum = 0.0;
    for ( size_t ititMoved(0), itFixed(0);
          itFixed < this->_FixedIndices->length();
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
    LOG("Number of moveable indices = {}" , this->_MoveableIndices->length() );
    ASSERT_gt(this->_MoveableIndices->length(),0);
    LOG("sumOfSquares = {}" , sumOfSquares );
    double meanSumOfSquares = sumOfSquares/this->_MoveableIndices->length();
    LOG("meanSumOfSquares = {}" , meanSumOfSquares );
    double rms = sqrt(meanSumOfSquares);
    LOG("rms = {}" , rms );
    return rms;
}

CL_LISPIFY_NAME("setFixedPoints");
CL_DEFMETHOD void SuperposeEngine_O::setFixedPoints(core::ComplexVector_byte32_t_sp fi, geom::SimpleVectorCoordinate_sp fc) {
    LOG("SuperposeEngine_O::setFixedPoints --> number of points={}", fc->length());
    if (!((fi->length() >= 3) && ((fi->length() <= fc->length())))) {
      SIMPLE_ERROR("In SuperposeEngine setFixedPoints - number of indices fi->length() = {} must be at least 3\nNumber of "
                   "coordinates --> fc->length() = {} must be >= to fi->size({})",
                   fi->length(), fc->length(), fi->length());
    }
    this->_FixedIndices = gc::As<core::ComplexVector_byte32_t_sp>(fi);
    geom::SimpleVectorCoordinate_sp data = fc->copy();
    this->_FixedCoordinates = geom::ComplexVectorCoordinate_O::make_vector(data);
}

CL_LISPIFY_NAME("setFixedAllPoints");
CL_DEFMETHOD void        SuperposeEngine_O::setFixedAllPoints( core::Array_sp afc )
{
  if (gc::IsA<geom::SimpleVectorCoordinate_sp>(afc) ) {
    auto fc = gc::As_unsafe<geom::SimpleVectorCoordinate_sp>(afc);
    this->_FixedCoordinates = geom::ComplexVectorCoordinate_O::make_vector(fc->length(),Vector3(),core::make_fixnum(fc->length()));
    this->_FixedIndices = core::ComplexVector_byte32_t_O::make_vector(fc->length());
    size_t ii;
    for ( ii=0; ii<fc->length(); ii++ ) {
      (*this->_FixedCoordinates)[ii] = (*fc)[ii];
      (*this->_FixedIndices)[ii] = ii;
//    printf("%s:%d Set fixed coordinate[%zu] -> %lf, %lf, %lf\n", __FILE__, __LINE__, ii, (*fc)[ii].getX(), (*fc)[ii].getY(), (*fc)[ii].getZ());
      LOG("setFixedAllPoints index@{} " , ii  );
    }
    return;
  } else if (gc::IsA<NVector_sp>(afc) ) {
    auto fc = gc::As_unsafe<NVector_sp>(afc);
    ASSERTF(fc->length()>=3*3,("You must have at least three moveable points and there are only {}") , fc->length()/3 );
    this->_FixedCoordinates = geom::ComplexVectorCoordinate_O::make_vector(fc->length()/3,Vector3(),core::make_fixnum(fc->length()/3));
    this->_FixedIndices = core::ComplexVector_byte32_t_O::make_vector(fc->length()/3);
    uint                ii;
    for ( ii=0; ii<fc->length()/3; ii++ )
    {
      (*this->_FixedCoordinates)[ii].set((*fc)[(ii*3)+0], (*fc)[(ii*3)+1], (*fc)[(ii*3)+2] );
      (*this->_FixedIndices)[ii] = ii;
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
    ASSERTP( (mi->length()>=3), "There must be at least three indices" );
    ASSERTP( mi->length()<=mc->length(), "There must be at least as many coordinates as indices");
    this->_MoveableIndices = gc::As<core::ComplexVector_byte32_t_sp>(mi);
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
    this->_MoveableIndices = core::ComplexVector_byte32_t_O::make_vector(mc->length());
    for ( ii=0; ii<mc->length(); ii++ )
    {
      (*this->_MoveableCoordinates)[ii] = (*mc)[ii];
      (*this->_MoveableIndices)[ii] = ii;
      LOG("setMoveableAllPoints index@{} " , ii  );
//    printf("%s:%d Set moveable coordinate[%d] -> %lf, %lf, %lf\n", __FILE__, __LINE__, ii, (*mc)[ii].getX(), (*mc)[ii].getY(), (*mc)[ii].getZ());
    }
    return;
  } else if (gc::IsA<NVector_sp>(amc)) {
    auto mc = gc::As_unsafe<NVector_sp>(amc);
    ASSERTF(mc->length()>=3*3,("You must have at least three moveable points and there are only {}") , mc->length() );
    this->_MoveableCoordinates = geom::ComplexVectorCoordinate_O::make_vector(mc->length()/3,Vector3(),core::make_fixnum(mc->length()));
    this->_MoveableIndices = core::ComplexVector_byte32_t_O::make_vector(mc->length()/3);
    uint                ii;
    for ( ii=0; ii<mc->length()/3; ii++ )
    {
      (*this->_MoveableCoordinates)[ii].set((*mc)[(ii*3)+0], (*mc)[(ii*3)+1], (*mc)[(ii*3)+2] );
      (*this->_MoveableIndices)[ii] = ii;
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
CL_DEFUN void chem__superpose_one(Matter_sp moveable_matter, core::List_sp moveable_atoms, core::List_sp fixed_atoms)
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
  geom::SimpleVectorCoordinate_sp coords_fixed = chem__make_simple_vector_coordinate_from_atom_list(fixed_atoms);
  geom::SimpleVectorCoordinate_sp coords_moveable = chem__make_simple_vector_coordinate_from_atom_list(moveable_atoms);
  SuperposeEngine_sp engine = SuperposeEngine_O::create();
  engine->setFixedAllPoints(coords_fixed);
  engine->setMoveableAllPoints(coords_moveable);
  Matrix transform = engine->superpose();
  moveable_matter->applyTransformToAtoms(transform);
}

};
