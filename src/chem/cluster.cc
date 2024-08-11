/*
    File: cluster.cc
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

#include <cando/chem/zMatrix.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/spanningLoop.h>
#include <cando/chem/loop.h>
#include <clasp/core/wrappers.h>
#include <cando/chem/cluster.h>

#include<iostream> 
#include<list>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>


#ifdef USING_OMP
#include<omp.h>
#endif

namespace chem {
Kmeans_O::Kmeans_O(int k)
{
	_K = k;

}
Kmeans_O::Kmeans_O(int k,int pointnumber)
{
  _K=k;
  _PointNumber = pointnumber;
}

// this function is used for compareing calculation time (Using omp)
void Kmeans_O::InitSpecifiedCenters(core::SimpleVector_sp centers)
{
  for (int i = 0; i < _K; i++)
  {
    (*centers)[i] = (*this->_Points)[i];
  }

}
CL_DEFMETHOD void  Kmeans_O::InitCenters(core::SimpleVector_sp centers)
{
#ifdef USING_OMP
#pragma ompparallel for 
#endif
  std::map<int, int> uniqueMap;
  while (uniqueMap.size() < _K)
  {
    int id = my_thread->random() % this->_PointNumber;
    uniqueMap.insert(std::pair<int, int>(id, id));
  }
  std::map<int, int>::iterator itMap = uniqueMap.begin();
  for (int i = 0; i < _K; i++)
  {

    int id = itMap->first;
    size_t it = 0;
    int count = 0;
    while (count != id)
    {
      it++;
      count++;
    }
    (*centers)[i] = (*this->_Points)[it];
    itMap++;
  }
}

/* Add the points to the Kmeans object
 */
void Kmeans_O::InitPoints(core::SimpleVector_sp points)
{
  this->_Points = points;
  this->_PointNumber = core::cl__length(points);
}

/* Return the index of the nearest center to the point
 */
int   Kmeans_O::NearestCenter(core::SimpleVector_sp centers, Point p)
{
  float minDistance = std::numeric_limits<float>::max();
  int k_id = -1;
  float dis;
  if (_K==0) SIMPLE_ERROR("There are no centers to compare to point _K == {}", _K );
  for (int k = 0; k < _K; k++ )
  {
    dis = Distance<floatType>(p, gc::As<Point>((*centers)[k]) );
    ENSURE_NOT_NAN(dis);
    if (dis < minDistance) {
      minDistance = dis;
      k_id = k;
    }
  }
  return k_id;
}

CL_DOCSTRING("Return the point at index idx");
CL_DEFMETHOD core::T_sp Kmeans_O::GetPoint(size_t idx) {
  if (idx<this->_Points->length()) {
    return (*(this->_Points))[idx];
  }
  SIMPLE_ERROR("Index {} out of range of number of points {}d", idx, this->_Points->length());
}

CL_DOCSTRING("Return all the points");
CL_DEFMETHOD core::SimpleVector_sp Kmeans_O::GetPoints() const {
  return this->_Points;
}


/* Assign every point to a cluster
*/
void  Kmeans_O::Cluster(core::SimpleVector_sp centers, Clusters clusters)
{
//  core::List_sp pointsIter = this->_Points;
  size_t pointsIter(0);
  for (int p = 0; p < _PointNumber; p++, pointsIter++)
  {
///    pointsIter = CONS_CDR(pointsIter);
    int nearest = NearestCenter(centers, gc::As<Point>((*this->_Points)[pointsIter]));
    if (nearest<0 || clusters->length()<=nearest) SIMPLE_ERROR("NearestCenter {} returned an index that is out of bounds ({})",
                                                              nearest, clusters->length());
    (*clusters)[p] = NearestCenter(centers, gc::As<Point>((*this->_Points)[pointsIter]));
  }
}

/* Calculate the center of each cluster.
 * Return the number of empty clusters.
*/
__attribute__((optnone))
CL_DEFMETHOD int Kmeans_O::Center(core::SimpleVector_sp centers, Clusters clusters)
{
  int emptyClusters = 0;
  core::T_sp one = (*this->_Points)[0];
  size_t point_length = gc::As<Point>(one)->length();
  for (int kk = 0; kk < this->_K; ++kk ) {
    (*centers)[kk] = core::SimpleVector_float_O::make(point_length,0.0,true);
  }
  std::vector<int> count(this->_K, 0);
  size_t centerIter(0);
  size_t pointsIter(0);
//  core::List_sp pointsIter = this->_Points;
//  core::List_sp centerIter = this->_Centers;
  for (int p = 0; p < _PointNumber; p++, pointsIter++ )
  {
//    pointsIter = CONS_CDR(pointsIter);
    size_t clusterIndex = (*clusters)[p];
    size_t centersSize = centers->length();
    if (centersSize<=clusterIndex) SIMPLE_ERROR("clusterIndex@{} is {}, which is out of bounds ({})", (void*)&((*clusters)[p]), clusterIndex, centersSize );
    core::T_sp tdest = (*centers)[clusterIndex];
    core::T_sp tsource = (*this->_Points)[pointsIter];
    if (!gc::IsA<Point>(tdest)) SIMPLE_ERROR("Dest is not a Point object it is {}", (void*)tdest.raw_());
    if (!gc::IsA<Point>(tsource)) SIMPLE_ERROR("Source is not a Point object it is {}", (void*)tsource.raw_());
    Point dest = gc::As<Point>(tdest);
    Point source = gc::As<Point>(tsource);
    AddToPoint(dest,source);
    count[(*clusters)[p]]++;
  }
  for (int i = 0; i < centers->length();i++)
  {
    if (count[i]>0) {
      ScalePoint(gc::As<Point>((*centers)[i]),1.0/count[i]);
    } else {
      (*centers)[i] = nil<core::T_O>(); // Set center to NIL
      emptyClusters++;
    }
    centerIter++;
  }
  return emptyClusters;
}
/*
 * Return a vector of the closest points in each cluster to the cluster center
 */
CL_DEFMETHOD core::SimpleVector_sp Kmeans_O::CenterOnPoints(core::SimpleVector_sp centers, Clusters clusters)
{
  core::SimpleVector_sp centerPoints = this->EmptyCenters();
  std::vector<float> closestDistances(this->_K,std::numeric_limits<float>::max());
  for (int pi=0; pi<this->_Points->length(); pi++ ) {
    int cluster = (*clusters)[pi];
    float distToCenter = Distance<floatType>(gc::As<Point>((*centers)[cluster]),gc::As<Point>((*this->_Points)[pi]));
    if (distToCenter < closestDistances[cluster]) {
      closestDistances[cluster] = distToCenter;
      (*centerPoints)[cluster] = (*this->_Points)[pi];
    }
  }
  return centerPoints;
}

/*
 * Return a vector of the closest points in each cluster to the cluster center
 */
CL_DEFMETHOD core::SimpleVector_float_sp Kmeans_O::PointDistancesToClusterCenters(core::SimpleVector_sp centers, Clusters clusters)
{
  core::SimpleVector_float_sp pointDistances = core::SimpleVector_float_O::make(this->_Points->length(),std::numeric_limits<float>::max(),true);
  for (int pi=0; pi<this->_Points->length(); pi++ ) {
    int cluster = (*clusters)[pi];
    float distToCenter = Distance<floatType>(gc::As<Point>((*centers)[cluster]),gc::As<Point>((*this->_Points)[pi]));
    (*pointDistances)[pi] = distToCenter;
  }
  return pointDistances;
}

/* Return a SimpleVector_int32_t_sp that store the cluster index for each point.
 */
CL_DEFMETHOD Clusters Kmeans_O::EmptyClusters() {
  Clusters clusters = core::SimpleVector_int32_t_O::make(this->_Points->length());
  return clusters;
}

CL_DEFMETHOD core::SimpleVector_sp Kmeans_O::EmptyCenters() {
  core::SimpleVector_sp centers = core::SimpleVector_O::make(this->_K);
  return centers;
}

/* Run the Kmean algorithm.
 * Return the number of times the algorithm was restarted.
 */
CL_DEFMETHOD int Kmeans_O::RunKmean(core::SimpleVector_sp centers, Clusters clusters, bool centerOnPoints )
{
  if (centers->length() != this->_K) {
    SIMPLE_ERROR("There is a mismatch between the number of centers {} and _K {}", centers->length(), this->_K);
  }
  if (clusters->length() != this->_Points->length()) {
    SIMPLE_ERROR("There is a mismatch between the number of clusters {} and number of points {}", clusters->length(), this->_Points->length());
  }
  core::SimpleVector_sp originalCenters = centers;
  core::SimpleVector_sp oldCenter;
  this->_MaxIteration = 100;
  int restarts = 0;
 top:
  for (int iteration = 0; iteration < this->_MaxIteration; iteration++)
  {
    // printf("%s:%d:%s iteration = %d\n", __FILE__, __LINE__, __FUNCTION__, iteration );
    oldCenter = centers;
    //printf("%s:%d:%s oldCenter = %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(oldCenter).c_str());
		//PrintPointLis(oldCenter);
    this->Cluster(centers,clusters);
    //printf("%s:%d:%s clusters= %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(clusters).c_str());
    centers = this->EmptyCenters();
    int emptyClusters = this->Center(centers,clusters);
    if (emptyClusters>0) {
      this->InitCenters(centers); // Initialize the centers and try again.
      restarts++;
      if (restarts>2) goto fail;
      goto top;
    }
    if (centerOnPoints) centers = this->CenterOnPoints(centers,clusters);
		//std::cout << "-------------------------\n";
    //printf("%s:%d:%s _Center = %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(this->_Centers).c_str());
		//PrintPointLis(_Centers);
    float sum = 0;
    // std::list<Point>::iterator currentIter = _Centers.begin();
    // std::list<Point>::iterator oldIter = oldCenter.begin();
    size_t currentIter = 0;
    size_t oldIter = 0;
    for (int k = 0; k < _K;k++)
    {
      sum += Distance<floatType>(gc::As<Point>((*oldCenter)[oldIter]), gc::As<Point>((*centers)[currentIter]));
      oldIter++;
      currentIter++;
    }
    //std::cout << "iteration "<< iteration<<" sum " << sum << std::endl;
    if (sum < 0.0001)
    {
      break;
    }
  }
  for ( int k=0; k< this->_K; k++ ) {
    (*originalCenters)[k] = (*centers)[k];
  }
  return restarts;
 fail:
  return restarts;
}
 
KmeansPlusPlus_O::KmeansPlusPlus_O(int k, int pointnumber) :Kmeans_O(k, pointnumber)
{
	
}
KmeansPlusPlus_O::KmeansPlusPlus_O(int k) : Kmeans_O(k)
{

}

int KmeansPlusPlus_O::NearestCenterDistance2(core::SimpleVector_sp centers, Point p, int numCenters, float &minDistance2)
{
  minDistance2 = std::numeric_limits<float>::max();
  int k_id = -1;
  float dis2;
	// std::list<Point>::iterator centersIter = _Centers.begin();
  size_t centersIter(0);
#ifdef USING_OMP
#pragma ompparallel for 
#endif
  for (int kk = 0; kk < numCenters; kk++ ) {
    core::T_sp tpnt = (*centers)[kk];
    if (tpnt.nilp()) SIMPLE_ERROR("The center at index %d is NIL", kk );
    Point pnt = gc::As<Point>(tpnt);
    dis2 = Distance2<floatType>(p, pnt );
    if (dis2 < minDistance2) {
      minDistance2 = dis2;
      k_id = kk;
    }
  }
  return k_id;
}


CL_DEFMETHOD void KmeansPlusPlus_O::InitCenters(core::SimpleVector_sp centers)
{
//  printf("%s:%d:%s  K = %d\n", __FILE__, __LINE__, __FUNCTION__, this->_K );
  size_t it = my_thread->random() % this->_PointNumber;
  (*centers)[0] = (*this->_Points)[it];
  std::vector<float> nearestDis2(_PointNumber,0.0);
  for ( int kk = 1; kk < _K; kk++ ) {
    float min_distance2;
    float sum = 0.0;
    std::vector<float>::iterator dis2It = nearestDis2.begin();
    for (int p = 0; p < _PointNumber; p++) {
      core::T_sp tpnt = (*this->_Points)[p];
      Point pnt = gc::As<Point>(tpnt);
      int nearest = NearestCenterDistance2(centers, pnt, kk, min_distance2);
      *dis2It = min_distance2;
//      printf("%s:%d:%s p=%d/%d *dist2It = %f sum= %f\n", __FILE__, __LINE__, __FUNCTION__, p, _PointNumber, *dis2It, sum );
      sum += *dis2It;
      dis2It++;
    }
//    printf("%s:%d:%s final sum = %f\n", __FILE__, __LINE__, __FUNCTION__, sum );
    dis2It = nearestDis2.begin();
    size_t it = my_thread->random() & 0xFFFFFFFF;
    float div = (float)0x100000000;
    float frac = it/div;
    float probability = sum*frac;
//    printf("%s:%d:%s  it = %lu frac = %f probability = %f\n", __FILE__, __LINE__, __FUNCTION__, it, frac, probability );
    sum = 0.0;
    for (int p=0; p < _PointNumber; p++ ) {
      sum += *dis2It;
//      printf("%s:%d:%s  sum= %f probability = %f\n", __FILE__, __LINE__, __FUNCTION__, sum, probability );
      if ( probability < sum ) {
//        printf("%s:%d:%s  match sum= %f\n", __FILE__, __LINE__, __FUNCTION__, sum );
        (*centers)[kk] = (*this->_Points)[p];
        break;
      }
      dis2It++;
    }
  }
}

CL_DEFMETHOD float Kmeans_O::Wcss(core::SimpleVector_sp centers, Clusters clusters) {
  std::vector<float> sums2(centers->length(),0.0);
  for ( int ii=0; ii < this->_Points->length(); ii++ ) {
    int cluster = (*clusters)[ii];
    sums2[cluster] += Distance2<floatType>(gc::As<Point>((*this->_Points)[ii]),gc::As<Point>((*centers)[cluster]));
  }
  return std::accumulate(sums2.begin(),sums2.end(),0.0);
}

CL_DEFMETHOD core::T_mv Kmeans_O::SilhouetteCoefficient(Clusters clusters) {
  int k = this->_K;
  Point ss = core::SimpleVector_float_O::make(this->_Points->length(),0.0,true);
  vector<int> counts(k, 0);
  for (int i = 0; i< this->_Points->length(); i++ ) {
    counts[(*clusters)[i]]++;
  }
  auto centers = this->EmptyCenters();
  this->Center(centers,clusters);
  for (int i = 0; i < this->_Points->length(); i++) {
    double a = 0.0, b = std::numeric_limits<double>::max();
    for (int j = 0; j < this->_Points->length(); j++) {
      if (i == j) continue;
      if ((*clusters)[i] == (*clusters)[j]) {
        float d = Distance<floatType>(gc::As<Point>((*this->_Points)[i]),gc::As<Point>((*this->_Points)[j]));
        a += d;
      } else {
        b = fmin(b, Distance<floatType>(gc::As<Point>((*this->_Points)[i]),gc::As<Point>((*this->_Points)[j])));
      }
    }
    if (counts[(*clusters)[i]] > 1) {
      a /= (counts[(*clusters)[i]] - 1);
    }
    (*ss)[i] = (b - a) / fmax(a, b);
  }
  double sum_s = 0.0;
  for (int i = 0; i < this->_Points->length(); i++) {
    sum_s += (*ss)[i];
  }
  return Values(core::make_single_float(sum_s/this->_Points->length()), ss );
}




CL_DEFUN KmeansPlusPlus_sp chem__make_kmeans_PLUS__PLUS_(size_t k, core::SimpleVector_sp points) {
  auto kpp = gctools::GC<KmeansPlusPlus_O>::allocate(k);
  kpp->InitPoints(points);
  return kpp;
};


};


namespace chem {

// Inspired by https://github.com/james-yoo/DBSCAN/blob/master/dbscan.cpp


int DBSCAN::run()
{
  int clusterID = 1;
  for(auto iter = 0; iter<this->_Points->length(); iter++ ) {
    if ( (*this->_Clusters)[iter] == UNCLASSIFIED ) {
      if ( this->expandCluster( iter, clusterID) != FAILURE ) {
        clusterID += 1;
      }
    }
  }
  return 0;
}

int DBSCAN::expandCluster(size_t point_index, int clusterId )
{
  Point point = gc::As<Point>((*this->_Points)[point_index]);
  std::vector<int> clusterSeeds = this->calculateCluster(point_index);
  if (clusterSeeds.size()<this->_minPoints) {
    (*this->_Clusters)[point_index] = NOISE;
    return FAILURE;
  } else {
    int index = 0;
    int indexCorePoint = 0;
    for (size_t iterSeeds = 0; iterSeeds<clusterSeeds.size(); iterSeeds++ ) {
      (*this->_Clusters)[iterSeeds] = clusterId;
      if (point_index == iterSeeds) {
        indexCorePoint = index;
      }
      ++index;
    }
    clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint) ;
    for( vector<int>::size_type i = 0, n = clusterSeeds.size() ; i < n; ++i ) {
      vector<int> clusterNeighors = this->calculateCluster(clusterSeeds[i]);
      if ( clusterNeighors.size() >= this->_minPoints ) {
        vector<int>::iterator iterNeighbors ;
        for ( iterNeighbors = clusterNeighors.begin() ; iterNeighbors != clusterNeighors.end(); ++iterNeighbors ) {
          if ((*this->_Clusters)[*iterNeighbors] == UNCLASSIFIED || (*this->_Clusters)[*iterNeighbors] == NOISE ) {
            if ((*this->_Clusters)[*iterNeighbors] == UNCLASSIFIED ) {
              clusterSeeds.push_back(*iterNeighbors) ;
              n = clusterSeeds.size();
            }
            (*this->_Clusters)[*iterNeighbors] = clusterId;
          }
        }
      }
    }
    return SUCCESS;
  }
}


std::vector<int> DBSCAN::calculateCluster(size_t point_index)
{
  Point point = gc::As<Point>((*this->_Points)[point_index]);
  int index = 0;
  std::vector<int> clusterIndex;
  for ( size_t iter = 0; iter< this->_Points->length(); iter++ ) {
    Point tpoint = gc::As<Point>((*this->_Points)[iter]);
    if (Distance(point,tpoint)<=this->_epsilon) {
      clusterIndex.push_back(iter);
    }
  }
  return clusterIndex;
}

#if 0
CL_DEFUN DBScan_sp chem__make_DBScan(core::SimpleVector_sp points) {
  auto kpp = gctools::GC<DBScan_O>::allocate(points);
  // Initialize DBScan object
  return kpp;
};
#endif

};
