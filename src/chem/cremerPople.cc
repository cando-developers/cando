/*
 *      File:   cremerPople.c
 *
 * Calculate Cremer-Pople parameters for rings and convert Cremer-Pople parameters into ring coordinates.
 *
 * See:
 *      General definition of ring puckering coordinates
 *      D. Cremer and J. A. Pople
 *      Journal of the American Chemical Society 1975 97 (6), 1354-1358
 *      DOI: 10.1021/ja00839a011
 *
 *      Understanding Ring Puckering in Small Molecules and Cyclic Peptides
 *      Lucian Chan, Geoffrey R. Hutchison, and Garrett M. Morris
 *      Journal of Chemical Information and Modeling 2021 61 (2), 743-755
 *      DOI: 10.1021/acs.jcim.0c01144
 */

#include <stdio.h>
#include <clasp/core/foundation.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/loop.h>
#include <cando/geom/ovector3.h>
#include <cando/geom/omatrix.h>
#include <clasp/core/bformat.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/aggregate.h>



namespace chem {

CL_DEFUN
string chem__cpring_mathematica(NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3) {
  size_t N = indexes3->length();
  stringstream ss;
  ss << fmt::format("{{ ");
  for ( size_t j=0; j<N-1; j++ ) {
    size_t index3 = (*indexes3)[j];
    ss << fmt::format("{{ {}, {}, {} }},\n",
                      (*coordinates)[index3+0],
                      (*coordinates)[index3+1],
                      (*coordinates)[index3+2] );
  }
  size_t index3 = (*indexes3)[N-1];
  ss << fmt::format("{{ {}, {}, {} }} }}\n",
                    (*coordinates)[index3+0],
                    (*coordinates)[index3+1],
                    (*coordinates)[index3+2] );
  return ss.str();
};

CL_DEFUN
core::T_mv chem__cpring_5_sample() {
  NVector_sp vec;
  vec = NVector_O::create(15);
  (*vec)[ 0] = 0.00000000;
  (*vec)[ 1] = 0.00000000;
  (*vec)[ 2] = 0.00000000;
  (*vec)[ 3] = 1.39513100;
  (*vec)[ 4] = 0.39017000;
  (*vec)[ 5] = 0.50319800;
  (*vec)[ 6] = 2.24910500;
  (*vec)[ 7] = -0.83685300;
  (*vec)[ 8] = 0.15180300;
  (*vec)[ 9] = 1.29665900;
  (*vec)[10] = -2.04400400;
  (*vec)[11] = 0.24315500;
  (*vec)[12] = -0.11423700;
  (*vec)[13] = -1.46402700;
  (*vec)[14] = 0.44433700;
  auto indexes3 = core::SimpleVector_byte32_t_O::make(5);
  (*indexes3)[0] = 0;
  (*indexes3)[1] = 3;
  (*indexes3)[2] = 6;
  (*indexes3)[3] = 9;
  (*indexes3)[4] = 12;
  return Values( vec, indexes3 );
}


CL_DEFUN
core::T_mv chem__cpring_6_sample() {
  NVector_sp vec;
  vec = NVector_O::create(18);
  (*vec)[ 0] =    0.00000000;
  (*vec)[ 1] =    0.00000000;
  (*vec)[ 2] =    0.00000000;
  (*vec)[ 3] =    0.02975530;
  (*vec)[ 4] =    -1.46190471;
  (*vec)[ 5] =    -0.46914502;
  (*vec)[ 6] =    1.30976685;
  (*vec)[ 7] =    -2.16378982;
  (*vec)[ 8] =    0.00737943;
  (*vec)[ 9] =    2.55854520;
  (*vec)[10] =    -1.40952475;
  (*vec)[11] =    -0.47197914;
  (*vec)[12] =    2.52878990;
  (*vec)[13] =    0.05237997;
  (*vec)[14] =    -0.00283412;
  (*vec)[15] =    1.24877835;
  (*vec)[16] =    0.75426508;
  (*vec)[17] =    -0.47935857;
  auto indexes3 = core::SimpleVector_byte32_t_O::make(6);
  (*indexes3)[0] = 0;
  (*indexes3)[1] = 3;
  (*indexes3)[2] = 6;
  (*indexes3)[3] = 9;
  (*indexes3)[4] = 12;
  (*indexes3)[5] = 15;
  return Values( vec, indexes3 );
}



CL_DEFUN
Vector3 chem__cpring_center_of_geometry( NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3 ) {
  Vector3 sum(0.0,0.0,0.0);
  size_t count = 0;
  for ( size_t idx=0; idx<indexes3->length(); idx++ ) {
    size_t index3 = (*indexes3)[idx];
    Vector3 pos(coordinates,index3,Unsafe());
    sum = sum + pos;
    count++;
  }
  Vector3 center = sum.multiplyByScalar(1.0/count);
  return center;
}

#if 0
CL_DEFUN
NVector_sp chem__cpring_translate_to_center_of_geometry( NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3 ) {
  if (coordinates->length()==0) {
    SIMPLE_ERROR("The coordinates must not be length=0");
  }
  if (coordinates->length()%3 != 0) {
    SIMPLE_ERROR("The coordinates must have a multiple of 3 elements");
  }
  Vector3 center = chem__cpring_center_of_geometry(coordinates,indexes3);
  NVector_sp trans = NVector_O::create(coordinates->length());
  for ( size_t idx=0; idx<coordinates->length(); idx+=3 ) {
    Vector3 pos(coordinates,idx,Unsafe());
    pos = pos - center;
    (*trans)[idx]   = pos.getX();
    (*trans)[idx+1] = pos.getY();
    (*trans)[idx+2] = pos.getZ();
  }
  return trans;
}
#endif

NVector_sp coordinates_to_zjs( NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3 ) {
  Vector3 center = chem__cpring_center_of_geometry( coordinates, indexes3 ); 
  size_t N = indexes3->length();
  Vector3 rp(0.0,0.0,0.0);
  Vector3 rpp(0.0,0.0,0.0);
  for (size_t j = 0; j<N; j++ ) {
    size_t idx = (*indexes3)[j];
    Vector3 rj(coordinates,idx,Unsafe());
    rj = rj - center; // translate to center of geometry
    double term = 2.0*M_PI*j/N;
    rp = rp + rj.multiplyByScalar(sin(term));
    rpp = rpp + rj.multiplyByScalar(cos(term));
  }
  //core::lisp_write(fmt::format("rp = {},{},{}\n",rp.getX(),rp.getY(),rp.getZ()));
  //core::lisp_write(fmt::format("rpp = {},{},{}\n",rpp.getX(),rpp.getY(),rpp.getZ()));
  Vector3 rpcrossrpp = rp.crossProduct(rpp);
  double len = rpcrossrpp.length();
  Vector3 n = rpcrossrpp.multiplyByScalar(1.0/len);
  //core::lisp_write(fmt::format("n = {},{},{}\n",n.getX(),n.getY(),n.getZ()));
  NVector_sp zjs = NVector_O::create(N);
  for (size_t j = 0; j<N; j++ ) {
    size_t idx = (*indexes3)[j];
    Vector3 rj(coordinates, idx, Unsafe());
    rj = rj - center;
    double zj = n.dotProduct(rj);
    (*zjs)[j] = zj;
  }
  return zjs;
};


CL_DEFUN
NVector_sp chem__cpring_coordinates_to_zjs( NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3 ) {
  NVector_sp zjs = coordinates_to_zjs( coordinates, indexes3 );
  return zjs;
}

NVector_sp ring_5_zjs_to_cremer_pople(NVector_sp zjs, double& phi2, double& q2) {
  size_t N = zjs->length();
  double cosSum = 0.0;
  double sinSum = 0.0;
  size_t mm = 2;
  for ( size_t j=0; j<N; j++ ) {
    double zj = (*zjs)[j];
    double term = 2.0*M_PI*mm*j/N;
    cosSum += zj*cos(term);
    sinSum += zj*sin(term);
    //printf("%s:%d:%s  j = %lu N = %lu  j<N = %d 8<5 = %d\n", __FILE__, __LINE__, __FUNCTION__, j, N, j<N, 8<5 );
  }
  double coef = sqrt(2.0/N);
  double q2cosphi2 =  coef*cosSum;
  double q2sinphi2 = -coef*sinSum;
  double sinphi2_over_cosphi2 = q2sinphi2/q2cosphi2;
  phi2 = atan(sinphi2_over_cosphi2);
  q2 = q2cosphi2/cos(phi2);
  return zjs;
}

CL_DEFUN
core::T_mv chem__cpring_5_coordinates_to_cremer_pople( NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3 ) {
  double phi2;
  double q2;
  NVector_sp zjs = coordinates_to_zjs( coordinates, indexes3 );
  ring_5_zjs_to_cremer_pople( zjs, phi2, q2 );
  return Values( core::DoubleFloat_O::create(phi2), core::DoubleFloat_O::create(q2) );
}

CL_DEFUN
core::T_mv chem__cpring_5_zjs_to_cremer_pople( NVector_sp zjs ) {
  double phi2;
  double q2;
  ring_5_zjs_to_cremer_pople( zjs, phi2, q2 );
  return Values( core::DoubleFloat_O::create(phi2), core::DoubleFloat_O::create(q2) );
}

void ring_6_zjs_to_cremer_pople(NVector_sp zjs, double& phi2, double& q2, double& qn2) {
  // Use the ring_5 func to calculate phi2 and q2
  ring_5_zjs_to_cremer_pople( zjs, phi2, q2);
  double inner = 0.0;
  size_t N = zjs->length();
  for ( size_t j=0; j<N; j++ ) {
    inner += pow(-1.0,j)*(*zjs)[j];
  }
  qn2 = 1.0/sqrt(N)*inner;
}

CL_DEFUN
core::T_mv chem__cpring_6_coordinates_to_cremer_pople( NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3 ) {
  double phi2;
  double q2;
  double qn2;
  NVector_sp zjs = coordinates_to_zjs( coordinates, indexes3 );
  ring_6_zjs_to_cremer_pople( zjs, phi2, q2, qn2 );
  return Values( core::DoubleFloat_O::create(phi2), core::DoubleFloat_O::create(q2), core::DoubleFloat_O::create(qn2) );
}

CL_DEFUN
core::T_mv chem__cpring_6_zjs_to_cremer_pople( NVector_sp zjs ) {
  double phi2;
  double q2;
  double qn2;
  ring_6_zjs_to_cremer_pople( zjs, phi2, q2, qn2 );
  return Values( core::DoubleFloat_O::create(phi2), core::DoubleFloat_O::create(q2), core::DoubleFloat_O::create(qn2) );
}



NVector_sp ring_5_cremer_pople_to_zjs(double phi2, double q2) {
  size_t nn = 5;
  NVector_sp zjs = NVector_O::create(nn);
  for ( size_t j = 0; j<nn; j++ ) {
    double term = phi2+2.0*M_PI*2.0*j/nn;
    double q2costerm = q2*cos(term);
    double zj = sqrt(2.0/nn)*q2costerm;
    (*zjs)[j] = zj;
  }
  return zjs;
}

NVector_sp ring_6_cremer_pople_to_zjs(double phi2, double q2, double qn2 ) {
  size_t nn = 6;
  NVector_sp zjs = NVector_O::create(nn);
  for ( size_t j = 0; j<nn; j++ ) {
    double term = phi2+2.0*M_PI*2.0*j/nn;
    double q2costerm = q2*cos(term);
    double zj = sqrt(2.0/nn)*q2costerm + 1.0/sqrt(nn)*qn2*pow(-1.0,j);
    (*zjs)[j] = zj;
  }
  return zjs;
}

CL_DEFUN
NVector_sp chem__cpring_5_cremer_pople_to_zjs( double phi2, double q2 ) {
  return ring_5_cremer_pople_to_zjs( phi2, q2 );
}


CL_DEFUN
NVector_sp chem__cpring_6_cremer_pople_to_zjs( double phi2, double q2, double qn2 ) {
  return ring_6_cremer_pople_to_zjs( phi2, q2, qn2 );
}

CL_DEFUN
void chem__cpring_coordinates_to_distances_and_angles( core::ComplexVector_double_sp distances, core::ComplexVector_double_sp angles, NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3 ) {
  size_t N = indexes3->length();
  for (size_t j = 0; j<N; j++ ) {
    int idxm1 = (j - 1);
    if (idxm1<0) idxm1 += N;
    int idx   = j;
    int idxp1 = j + 1;
    if (idxp1>=N) idxp1 -= N;
    size_t index3m1 = (*indexes3)[idxm1];
    size_t index3   = (*indexes3)[idx];
    size_t index3p1 = (*indexes3)[idxp1];
    Vector3 posm1((*coordinates)[index3m1],(*coordinates)[index3m1+1],(*coordinates)[index3m1+2]);
    Vector3 pos((*coordinates)[index3],(*coordinates)[index3+1],(*coordinates)[index3+2]);
    Vector3 posp1((*coordinates)[index3p1],(*coordinates)[index3p1+1],(*coordinates)[index3p1+2]);
    (*distances)[j] = geom::calculateDistance(pos,posp1);
    (*angles)[j] = geom::calculateAngle(posm1,pos,posp1);
  }
};


CL_DEFUN
core::T_mv chem__cpring_generate_coordinates( NVector_sp zjs, core::SimpleVector_double_sp r_values, core::SimpleVector_double_sp beta_rad_values, NVector_sp coordinates, core::SimpleVector_byte32_t_sp indexes3 ) {
  core::SimpleVector_double_sp r_prime_values = core::SimpleVector_double_O::create(zjs->length());
  core::SimpleVector_double_sp beta_prime_values = core::SimpleVector_double_O::create(zjs->length());
  for ( size_t ii = 0; ii<zjs->length(); ii++ ) {
    size_t jj = (ii + 1) % zjs->length();
    size_t kk = (ii + 2) % zjs->length();
    double rab = (*r_values)[ii];
    double za = (*zjs)[ii];
    double zb = (*zjs)[jj];
    double zc = (*zjs)[kk];
    double rabprime = fabs(sqrt(rab*rab-(za-zb)*(za-zb)));
    double rbc = (*r_values)[jj];
    double beta = (*beta_rad_values)[ii];
    double cosbeta_prime = (rab*rbc*cos(beta)-(za-zb)*(zc-zb))/(sqrt(rab*rab-(za-zb)*(za-zb))*sqrt(rbc*rbc-(zc-zb)*(zc-zb)));
    double beta_prime = acos(cosbeta_prime);
    (*r_prime_values)[ii] = rabprime;
    (*beta_prime_values)[ii] = beta_prime;
  }
  for ( size_t ii = 0; ii < zjs->length(); ii ++ ) {
    size_t index3 = (*indexes3)[ii];
    if (ii==0) {
      (*coordinates)[index3]   = 0.0;
      (*coordinates)[index3+1] = 0.0;
      (*coordinates)[index3+2] = (*zjs)[ii];
    } else if (ii==1) {
      size_t index3m1 = (*indexes3)[0];
      Vector3 posm1((*coordinates)[index3m1],(*coordinates)[index3m1+1],(*coordinates)[index3m1+2]);
      Vector3 pos = geom::geom__build_using_bond( (*r_prime_values)[ii-1], posm1 );
      (*coordinates)[index3]   = pos.getX();
      (*coordinates)[index3+1] = pos.getY();
      (*coordinates)[index3+2] = (*zjs)[ii];
    } else {
      size_t index3m2 = (*indexes3)[ii-2];
      size_t index3m1 = (*indexes3)[ii-1];
      Vector3 posm2((*coordinates)[index3m2+0],(*coordinates)[index3m2+1],(*coordinates)[index3m2+2]);
      Vector3 posm1((*coordinates)[index3m1+0],(*coordinates)[index3m1+1],(*coordinates)[index3m1+2]);
      Vector3 pos = geom::geom__build_using_bond_angle( (*r_prime_values)[ii-1], posm1, (*beta_prime_values)[ii-1], posm2 );
      (*coordinates)[index3]   = pos.getX();
      (*coordinates)[index3+1] = pos.getY();
      (*coordinates)[index3+2] = (*zjs)[ii];
    }
  }
  return Values( coordinates, r_prime_values, beta_prime_values );
}


};
