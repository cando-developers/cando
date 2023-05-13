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
  float sum = 0.0;
  for ( size_t ii=0; ii<dest->length(); ii++ ) {
    (*dest)[ii] += (*source)[ii];
  }
}

inline void ScalePoint(Point dest, float scale ) {
  for ( size_t ii=0; ii<dest->length(); ii++ ) {
    (*dest)[ii] *= scale;
  }
}


inline float Distance(Point p1, Point p2)
{
  float sum = 0.0;
  for ( size_t ii=0; ii<p1->length(); ii++ ) {
    float delta = (*p1)[ii] - (*p2)[ii];
    sum += delta * delta;
  };
  return std::sqrt(sum);
}


inline float Distance2(Point p1, Point p2)
{
  float sum = 0.0;
  for ( size_t ii=0; ii<p1->length(); ii++ ) {
    float delta = (*p1)[ii] - (*p2)[ii];
    sum += delta * delta;
  };
  return sum;
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
  int NearestCenter(core::SimpleVector_sp centers, Point p,int alreadyInitCenterNumber, float &minDistance);
private:

};

};



#endif

