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
#include <random>
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
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<int> distribution(0, _PointNumber - 1);
#ifdef USING_OMP
#pragma ompparallel for 
#endif
  std::map<int, int> uniqueMap;
  while (uniqueMap.size() < _K)
  {
    int id = distribution(gen);
    uniqueMap.insert(std::pair<int, int>(id, id));
  }
  std::map<int, int>::iterator itMap = uniqueMap.begin();
  for (int i = 0; i < _K; i++)
  {

    int id = itMap->first;
    size_t it;
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
  size_t centersIter(0);
  for (int k = 0; k < _K; k++, centersIter++ )
  {
    dis = Distance(p, gc::As<Point>((*centers)[k]) );
    if (dis < minDistance)
    {
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
  SIMPLE_ERROR("Index %ud out of range of number of points %ud", idx, this->_Points->length());
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
  float minDistance;
  for (int p = 0; p < _PointNumber; p++, pointsIter++)
  {
///    pointsIter = CONS_CDR(pointsIter);
    (*clusters)[p] = NearestCenter(centers, gc::As<Point>((*this->_Points)[pointsIter]));
  }
}

/* Calculate the center of each cluster.
 * Return the number of empty clusters.
*/
int  Kmeans_O::Center(core::SimpleVector_sp centers, Clusters clusters)
{
  int emptyClusters = 0;
  core::T_sp one = (*this->_Points)[0];
  size_t point_length = gc::As<Point>(one)->length();
  for (int kk = 0; kk < this->_K; ++kk ) {
    (*centers)[kk] = core::SimpleVector_float_O::make(point_length,0.0,true);
  }
  std::vector<int> count(this->_K, 0);
  size_t pointsIter(0);
  size_t centerIter(0);
//  core::List_sp pointsIter = this->_Points;
//  core::List_sp centerIter = this->_Centers;
  for (int p = 0; p < _PointNumber; p++, pointsIter++ )
  {
//    pointsIter = CONS_CDR(pointsIter);
    AddToPoint(gc::As<Point>((*centers)[(*clusters)[p]]),gc::As<Point>((*this->_Points)[pointsIter]));
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

core::SimpleVector_sp  Kmeans_O::ClosestPointIndexToCenter(Clusters clusters)
{
  core::SimpleVector_sp centers = this->EmptyCenters();
  this->Center(centers,clusters);
  core::SimpleVector_sp closestPointIndices = this->EmptyCenters();
  std::vector<float> closestDistances(this->_K,std::numeric_limits<float>::max());
  for (int pi=0; pi<this->_Points->length(); pi++ ) {
    int cluster = (*clusters)[pi];
    float distToCenter = Distance(gc::As<Point>((*centers)[cluster]),gc::As<Point>((*this->_Points)[pi]));
    if (distToCenter < closestDistances[cluster]) {
      closestDistances[cluster] = distToCenter;
      (*closestPointIndices)[cluster] = core::make_fixnum(pi);
    }
  }
  return closestPointIndices;
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
CL_DEFMETHOD int Kmeans_O::RunKmean(core::SimpleVector_sp centers, Clusters clusters)
{
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
      sum += Distance(gc::As<Point>((*oldCenter)[oldIter]), gc::As<Point>((*centers)[currentIter]));
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
int KmeansPlusPlus_O::NearestCenter(core::SimpleVector_sp centers, Point p, int alreadyInitCenterNumber, float &minDistance)
{
  minDistance = std::numeric_limits<float>::max();
  int k_id = -1;
  float dis;
	// std::list<Point>::iterator centersIter = _Centers.begin();
  size_t centersIter(0);
#ifdef USING_OMP
#pragma ompparallel for 
#endif
  for (int k = 0; k <= alreadyInitCenterNumber; centersIter++, k++)
  {
    dis = Distance(p, gc::As<Point>((*centers)[centersIter]));
    if (dis < minDistance)
    {
      minDistance = dis;
      k_id = k;
    }
  }
  return k_id;
}


CL_DEFMETHOD void KmeansPlusPlus_O::InitCenters(core::SimpleVector_sp centers)
{
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<int> distribution(0, _PointNumber - 1);
  int id = distribution(gen);
        // std::list<Point>::iterator it = _Points.begin();
  size_t it(id);
	// std::advance(it, id);
  for (int i = 0; i < _K; i++)
  {
    (*centers)[i] = (*this->_Points)[it];
  }

  float sum,min_distance;
 // std::list<Point>::iterator centersIter = _Centers.begin();
  // std::list<Point>::iterator pointIter = _Points.begin();
  size_t centersIter(0);
  size_t pointIter(0);
  std::list<float> nearestDis(_PointNumber,0);
  std::list<float>::iterator floatIt = nearestDis.begin();
	
  for (int k = 1; k < _K; centersIter++, k++)
  {
    sum = 0;
    pointIter = 0; //_Points.begin();
    floatIt = nearestDis.begin();
    for (int p = 0; p < _PointNumber; pointIter++, p++)
    {
      NearestCenter(centers, gc::As<Point>((*this->_Points)[pointIter]),k, min_distance);
      *floatIt = min_distance;
      sum += min_distance;
      floatIt++;
    }
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float probability = distribution(gen);
		//std::cout << "orignional sum " << sum << std::endl;
    sum = sum*probability;
		//std::cout << "orignional sum " << sum << std::endl;
    pointIter = 0; // _Points.begin();
    floatIt = nearestDis.begin();
    for (int p = 0; p < _PointNumber; pointIter++, floatIt++, p++)
    {
      sum =sum- *floatIt;
			 //std::cout <<"p " << p<< " sum "<<sum << std::endl;
      if (sum >0)
        continue;
      (*centers)[centersIter] = (*this->_Points)[pointIter];
      // centersIter->_x = pointIter->_x;
      // centersIter->_y = pointIter->_y;
      break;
    }
		 
  }
	/*float *d = malloc(sizeof(float)* len);

	point p, c;
	cent[0] = pts[rand() % len];
	for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
		sum = 0;
		for (j = 0, p = pts; j < len; j++, p++)
		{
			nearest(p, cent, n_cluster, d + j);
			sum += d[j];
		}
		sum = randf(sum);
		for (j = 0, p = pts; j < len; j++, p++)
		{
			if ((sum -= d[j]) > 0) continue;
			cent[n_cluster] = pts[j];
			break;
		}
	}
	for (j = 0, p = pts; j < len; j++, p++)
		p->group = nearest(p, cent, n_cluster, 0);*/
	
}
CL_DEFMETHOD float Kmeans_O::Wcss(core::SimpleVector_sp centers, Clusters clusters) {
  std::vector<float> sums2(centers->length(),0.0);
  for ( int ii; ii < this->_Points->length(); ii++ ) {
    int cluster = (*clusters)[ii];
    sums2[cluster] += Distance2(gc::As<Point>((*this->_Points)[ii]),gc::As<Point>((*centers)[cluster]));
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
        float d = Distance(gc::As<Point>((*this->_Points)[i]),gc::As<Point>((*this->_Points)[j]));
        a += d;
      } else {
        b = fmin(b, Distance(gc::As<Point>((*this->_Points)[i]),gc::As<Point>((*this->_Points)[j])));
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

