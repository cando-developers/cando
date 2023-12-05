#ifndef _KMEANS_H
#define _KMEANS_H

//
// Adapted from https://github.com/ieyjzhou/KmeansPlusPlus
// Originally written by  y.j.zhou.g@gmail.com.
// Made multidimensional by Christian Schafmeister
//


//#define  USING_OMP

#include<iostream> 
#include<list>
#include <chrono>
#include <random>

#include <clasp/core/common.h>
#include <clasp/core/ql.h>
#include <clasp/core/array_float.h>
#include <clasp/core/array_int32.h>


namespace chem {

FORWARD(Kmeans);
FORWARD(KmeansPlusPlus);

typedef core::SimpleVector_float_sp     Point;
typedef core::SimpleVector_int32_t_sp    Clusters;

inline void AddToPoint(Point dest, Point source ) {
  for ( size_t ii=0; ii<dest->length(); ii++ ) {
    (*dest)[ii] += (*source)[ii];
  }
}

inline void ScalePoint(Point dest, float scale ) {
  for ( size_t ii=0; ii<dest->length(); ii++ ) {
    (*dest)[ii] *= scale;
  }
}

template <typename Type=float>
inline float Distance2(Point p1, Point p2)
{
  float sum = 0.0;
  for ( size_t ii=0; ii<p1->length(); ii++ ) {
    float delta = (*p1)[ii] - (*p2)[ii];
    sum += delta * delta;
  };
  return sum;
}

template <typename Type=float>
inline float Distance(Point p1, Point p2)
{
  float sum2 = Distance2(p1,p2);
  return std::sqrt(sum2);
}

#define SIMD_CHUNK_SIZE 8
typedef float v8sf __attribute__((vector_size(SIMD_CHUNK_SIZE*sizeof(float))));
typedef float floatType; // can be float
//typedef v8sf floatType; // can be float


template <>
inline float Distance2<v8sf>(Point p1, Point p2)
{
  size_t chunks = p1->length()/SIMD_CHUNK_SIZE;
  size_t remainder = p1->length() - (chunks*SIMD_CHUNK_SIZE);
  float sum = 0.0;
  float* pf1 = &(*p1)[0];
  float* pf2 = &(*p2)[0];
  if (chunks) {
    for ( size_t chunk = 0; chunk<chunks; chunk+=SIMD_CHUNK_SIZE ) {
      // handle full chunks
      v8sf v1 = {pf1[chunk+0], pf1[chunk+1], pf1[chunk+2], pf1[chunk+3], pf1[chunk+4], pf1[chunk+5], pf1[chunk+6], pf1[chunk+7] };
      v8sf v2 = {pf2[chunk+0], pf2[chunk+1], pf2[chunk+2], pf2[chunk+3], pf2[chunk+4], pf2[chunk+5], pf2[chunk+6], pf2[chunk+7] };
      v8sf delta = v1 - v2;
      v8sf delta2 = delta*delta;
      sum += delta2[0]+ delta2[1]+ delta2[2]+ delta2[3]+ delta2[4]+ delta2[5]+ delta2[6]+ delta2[7];
    }
  }
  if (remainder) {
    v8sf vr1 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    v8sf vr2 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    size_t remainder_chunk = chunks*SIMD_CHUNK_SIZE;
    for ( size_t rem = 0; rem<remainder; rem++ ) {
      vr1[rem] = pf1[remainder_chunk+rem];
      vr2[rem] = pf2[remainder_chunk+rem];
    }
    v8sf delta = vr1 - vr2;
    v8sf delta2 = delta*delta;
    for ( size_t rem = 0; rem<remainder; rem++ ) {
      sum += delta2[rem];
    }
  }
  return sum;
}

template <>
inline float Distance<v8sf>(Point p1, Point p2)
{
  float sum2 = Distance2<v8sf>(p1,p2);
  return std::sqrt(sum2);
}


class Kmeans_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,Kmeans_O,"KMEANS",core::CxxObject_O);
public:
  Kmeans_O(int k, int pointnumber);
  Kmeans_O(int k);
  Kmeans_O() : _K(0) {};
  void InitPoints(core::SimpleVector_sp point_vector);
  virtual  void InitCenters(core::SimpleVector_sp centers);
  void InitSpecifiedCenters(core::SimpleVector_sp centers);
  int NearestCenter(core::SimpleVector_sp centers, Point p);
  Clusters EmptyClusters();
  core::SimpleVector_sp EmptyCenters();
  float Wcss(core::SimpleVector_sp centers, Clusters clusters );
  int RunKmean(core::SimpleVector_sp centers, Clusters clusters, bool centerOnPoints );
  void Cluster(core::SimpleVector_sp centers, Clusters clusters);
  int Center(core::SimpleVector_sp centers, Clusters clusters);
  core::SimpleVector_sp CenterOnPoints(core::SimpleVector_sp centers, Clusters clusters);
  core::SimpleVector_float_sp PointDistancesToClusterCenters(core::SimpleVector_sp centers, Clusters clusters);
  CL_DEFMETHOD int GetK() const {return this->_K;};
  core::SimpleVector_sp GetPoints() const;
  core::T_sp GetPoint(size_t idx);
  core::T_mv SilhouetteCoefficient(Clusters clusters);
public:
  core::SimpleVector_sp _Points;
  int _MaxIteration;
  int _K;
  int _PointNumber;
} ;

class KmeansPlusPlus_O: public Kmeans_O
{
  LISP_CLASS(chem,ChemPkg,KmeansPlusPlus_O,"KMEANS++",Kmeans_O);
public:
  KmeansPlusPlus_O(int k, int pointnumber);
  KmeansPlusPlus_O(int k);
  KmeansPlusPlus_O() : Kmeans_O(0) {};
  
  virtual void InitCenters(core::SimpleVector_sp centers);
  int NearestCenterDistance2(core::SimpleVector_sp centers, Point p,int alreadyInitCenterNumber, float &minDistance2);
private:

};


};



#endif

