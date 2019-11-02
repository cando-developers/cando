/*
    File: energyNonbond.cc
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
//#define	DEBUG_LEVEL_FULL

#include <clasp/core/foundation.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energyRigidBodyNonbond.h>
#include <cando/chem/rigidBodyEnergyFunction.h>
#include <clasp/core/array.h>
#include <clasp/core/evaluator.h>
#include <cando/geom/color.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/nativeVector.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/matter.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/nVector.h>
#include <cando/chem/loop.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

#if 0
#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
#define LOG_ENERGY BFORMAT_T
#endif

namespace chem
{

struct TypeParameters {
  double _Radius;
  double _Epsilon;
  TypeParameters(double radius, double epsilon) : _Radius(radius), _Epsilon(epsilon) {};
};

bool operator<(const TypeParameters& l, const TypeParameters& r) {
  return (l._Radius<r._Radius || (l._Radius==r._Radius && l._Epsilon<r._Epsilon));
}

CL_LISPIFY_NAME(energy-rigid-body-nonbond-number-of-nonbond-atoms);
CL_DEFMETHOD size_t EnergyRigidBodyNonbond_O::number_of_nonbond_atoms() const {
  return this->_AtomInfoTable.size();
}

CL_DEFMETHOD void EnergyRigidBodyNonbond_O::initializeCrossTerms(bool verbose)
{
  typedef std::map<TypeParameters,size_t> double_pair_map;
  double_pair_map types;
  std::vector<TypeParameters> typeParameters;
  size_t next_type_index = 0;
  for ( size_t i=0; i<this->_AtomInfoTable.size(); ++i ) {
    RigidBodyAtomInfo& info = this->_AtomInfoTable[i];
    TypeParameters key(info._Radius,info._Epsilon);
    double_pair_map::iterator it = types.find(key);
    if (it == types.end()) {
      types[key] = next_type_index;
      this->_AtomInfoTable[i]._TypeIndex = next_type_index;
      typeParameters.push_back(key);
      next_type_index++;
    } else {
      this->_AtomInfoTable[i]._TypeIndex = it->second;
    }
  }
  this->_NumberOfTypes = next_type_index;
  if (verbose) {
    core::write_bf_stream(BF("There are %d types\n") % this->_NumberOfTypes);
  }
  this->_CrossTerms.resize(next_type_index*next_type_index);
  double vdwScale = this->getVdwScale();
  for ( size_t xi=0; xi<next_type_index; xi++ ) {
    for (size_t yi=xi; yi<next_type_index; yi++ ) {
      TypeParameters& ea1 = typeParameters[xi];
      TypeParameters& ea2 = typeParameters[yi];
      double rStar = ea1._Radius+ea2._Radius;
      double epsilonij = sqrt(ea1._Epsilon*ea2._Epsilon);
      double rStar2 = rStar*rStar;
      double rStar6 = rStar2*rStar2*rStar2;
      double rStar12 = rStar6*rStar6;
      double dA = epsilonij*rStar12*vdwScale;
      double dC = 2.0*epsilonij*rStar6*vdwScale;
      RigidBodyNonbondCrossTerm crossTerm(dA,dC);
      if (verbose) {
        core::write_bf_stream(BF("Cross term for types %d - %d\n") % xi % yi );
        core::write_bf_stream(BF("    type: %d   radius-> %f   epsilon-> %f\n") % xi % ea1._Radius % ea1._Epsilon);
        core::write_bf_stream(BF("    type: %d   radius-> %f   epsilon-> %f\n") % yi % ea2._Radius % ea2._Epsilon);
        core::write_bf_stream(BF("    dA -> %f      dC -> %f\n") % dA % dC );
      }
      this->crossTerm(xi,yi) = crossTerm;
      this->crossTerm(yi,xi) = crossTerm;
    }
  }
}
          
EnergyRigidBodyNonbond_sp EnergyRigidBodyNonbond_O::make(core::Array_sp end_atoms) {
  if (end_atoms->length()<1) {
    SIMPLE_ERROR(BF("You must provide a vector of end atom indices with at least one end atom"));
  }
  core::SimpleVector_byte32_t_sp sv = core::SimpleVector_byte32_t_O::make(end_atoms->length());
  // Copy the atom indices
  size_t istart = 0;
  for ( size_t i(0), iEnd(end_atoms->length()); i<iEnd; ++i ) {
    Fixnum f = core::clasp_to_fixnum(end_atoms->rowMajorAref(i));
    if (istart>=f) {
      SIMPLE_ERROR(BF("The list of atom indices must be strictly increasing"));
    }
    (*sv)[i] = f;
  }
  GC_ALLOCATE_VARIADIC(EnergyRigidBodyNonbond_O,nb,sv);
  return nb;
}

void	EnergyRigidBodyNonbond_O::zeroEnergy()
{
  this->Base::zeroEnergy();
  this->_EnergyElectrostatic = 0.0;
  this->_EnergyVdw = 0.0;
}

double	EnergyRigidBodyNonbond_O::getEnergy()
{
  double	e;
  e = this->getVdwEnergy();
  e += this->getElectrostaticEnergy();
  return e;
}

 void EnergyRigidBodyNonbond_O::energyRigidBodyNonbondSetTerm(gc::Fixnum index, core::T_sp atom, double radius, double epsilon, double charge, const Vector3& position) {
  if (index < 0 || index >= this->_AtomInfoTable.size()) {
    SIMPLE_ERROR(BF("Index out of range %d - max is %d") % index % this->_AtomInfoTable.size());
  }
  if (this->_AtomInfoTable[index]._Object.boundp()) {
    SIMPLE_ERROR(BF("The rigid body nonbond term at %d has already been set to object %s") % index % _rep_(this->_AtomInfoTable[index]._Object));
  }
  RigidBodyAtomInfo info(atom,radius,epsilon,charge,position);
  this->_AtomInfoTable[index] = info;
}

CL_LISPIFY_NAME(energy-rigid-body-nonbond-get-position);
CL_DEFMETHOD
Vector3 EnergyRigidBodyNonbond_O::getPosition(size_t index)
{
  if (index < this->_AtomInfoTable.size()) {
    return this->_AtomInfoTable[index]._Position;
  }
  SIMPLE_ERROR(BF("index %lu is out of bounds") % index);
}


void	EnergyRigidBodyNonbond_O::dumpTerms()
{
  size_t istart = 0;
  for ( size_t rb = 0; rb<this->_RigidBodyEndAtom->length(); ++rb) {
    for ( size_t i = istart; i<(*this->_RigidBodyEndAtom)[rb]; ++i ) {
      RigidBodyAtomInfo& ai = this->_AtomInfoTable[i];
      core::write_bf_stream(BF("I1 = %3d  %s  %lf %lf %lf %s\n") % (rb*7) % _rep_(ai._Object) % ai._Radius % ai._Epsilon % ai._Charge % ai._Position.asString());
    }
    istart = (*this->_RigidBodyEndAtom)[rb];
  }
}

CL_DEFMETHOD core::List_sp EnergyRigidBodyNonbond_O::parts_as_list(NVector_sp pos)
{
  ql::list result;
  size_t istart = 0;
#undef	NONBOND_POSITION_RB_SET_PARAMETER
#define	NONBOND_POSITION_RB_SET_PARAMETER(x)	{}
#undef	NONBOND_POSITION_RB_SET_POSITION
#define	NONBOND_POSITION_RB_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_POSITION_RB_SET_POINT
#define	NONBOND_POSITION_RB_SET_POINT(x,ii,of)	{x=ii._Position.of;}
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termDeclares.cc>
#pragma clang diagnostic pop
  double am, bm, cm, dm, xm, ym, zm;
  double pxm, pym, pzm;
  int	I1;
  size_t I1start = 0;
  for ( size_t iI1 = 0; iI1<this->_RigidBodyEndAtom->length(); ++iI1 ) {
    ql::list helix;
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
    for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
      RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
      I1 = iI1*7;
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termCode.cc>
      helix << core::Cons_O::createList(core::DoubleFloat_O::create(plabmx),
                                        core::DoubleFloat_O::create(plabmy),
                                        core::DoubleFloat_O::create(plabmz),
                                        core::DoubleFloat_O::create(ea1._Radius));
    }
    I1start = I1end;
    result << helix.cons();
  }
  return result.cons();
}

struct coordinate_lookup {
  coordinate_lookup(NVector_sp coords) : _Coordinates(coords) {};
  NVector_sp _Coordinates;
  size_t _StartIndex;
  void setStartIndex(size_t startIndex) { this->_StartIndex = startIndex; };
  double getX() const { return (*this->_Coordinates)[this->_StartIndex+0];}
  double getY() const { return (*this->_Coordinates)[this->_StartIndex+1];}
  double getZ() const { return (*this->_Coordinates)[this->_StartIndex+2];}
};



CL_DOCSTRING(R"doc(This function is to test the rigid body transformation code.
Arguments:
rigid-body-pos : An nvector of N quaternion/origin coordinates (a,b,c,d,x,y,z).
end-indicesx3 : A vector of integer indices that index into 'coordinates' and indicate where
each object (defined by a contiguous set of points) in coordinates ends. The 'end-indicesx3'
values are indexes directly into 'coordinates', so they are x3.
coordinates : An nvector of coordinates.
output : A complex-vector-float where the transformed points are written.
)doc");
CL_DEFMETHOD core::ComplexVector_float_sp EnergyRigidBodyNonbond_O::write_rigid_body_coordinates_to_complex_vector_float(NVector_sp rigid_body_pos, core::Array_sp end_indicesx3, NVector_sp coordinates, core::ComplexVector_float_sp output)
{
  size_t istart = 0;
#undef	NONBOND_POSITION_RB_SET_PARAMETER
#define	NONBOND_POSITION_RB_SET_PARAMETER(x)	{}
#undef	NONBOND_POSITION_RB_SET_POSITION
#define	NONBOND_POSITION_RB_SET_POSITION(x,ii,of)	{x=rigid_body_pos->element(ii+of);}
#undef	NONBOND_POSITION_RB_SET_POINT
#define	NONBOND_POSITION_RB_SET_POINT(x,ii,of)	{x=ii.of;}
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termDeclares.cc>
#pragma clang diagnostic pop
  double am, bm, cm, dm, xm, ym, zm;
  double pxm, pym, pzm;
  int	I1;
  size_t coordIndex = 0;
  size_t I1start = 0;
  coordinate_lookup ea1(coordinates); 
  for ( size_t iI1 = 0; iI1<end_indicesx3->length(); ++iI1 ) {
    size_t I1end = core::clasp_to_size(end_indicesx3->rowMajorAref(iI1));
    for ( size_t I1cur = I1start; I1cur<I1end; I1cur+=3 ) {
      ea1.setStartIndex(I1cur);
      I1 = iI1*7;
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termCode.cc>
      output->vectorPushExtend((float)plabmx);
      output->vectorPushExtend((float)plabmy);
      output->vectorPushExtend((float)plabmz);
#if 0
      if (coordIndex==0) {
        core::write_bf_stream(BF("fill_pointer -> %d\n") % output->fillPointer() );
        core::write_bf_stream(BF("&plabmx -> %p\n") % (void*)&plabmx);
        core::write_bf_stream(BF("&output[0] -> %p\n") % (void*)&(*output)[0]);
        core::write_bf_stream(BF("Point[0] -> %f, %f, %f\n") % plabmx % plabmy % plabmz);
        core::write_bf_stream(BF("output[0] -> %f\n") % (*output)[0]);
        core::write_bf_stream(BF("hex dump output[0] -> %lX\n") % *(unsigned int*)&(*output)[0]);
        (*output)[0] = plabmx;
        core::write_bf_stream(BF("after output[0] -> %f\n") % (*output)[0]);
        core::write_bf_stream(BF("hex dump output[0] -> %lX\n") % *(unsigned int*)&(*output)[0]);
      }
      coordIndex += 3;
#endif
    }
    I1start = I1end;
  }
  return output;
}

CL_DEFMETHOD core::ComplexVector_float_sp EnergyRigidBodyNonbond_O::write_nonbond_atom_coordinates_to_complex_vector_float(NVector_sp pos, core::ComplexVector_float_sp parts)
{
  size_t istart = 0;
#undef	NONBOND_POSITION_RB_SET_PARAMETER
#define	NONBOND_POSITION_RB_SET_PARAMETER(x)	{}
#undef	NONBOND_POSITION_RB_SET_POSITION
#define	NONBOND_POSITION_RB_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_POSITION_RB_SET_POINT
#define	NONBOND_POSITION_RB_SET_POINT(x,ii,of)	{x=ii._Position.of;}
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termDeclares.cc>
#pragma clang diagnostic pop
  double am, bm, cm, dm, xm, ym, zm;
  double pxm, pym, pzm;
  int	I1;
  size_t I1start = 0;
  for ( size_t iI1 = 0; iI1<this->_RigidBodyEndAtom->length(); ++iI1 ) {
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
#if 0
    if (iI1<3) {
      core::write_bf_stream(BF("%s:%d I1start -> %d  I1end -> %d\n") % __FILE__ % __LINE__ % I1start % I1end);
    }
#endif
    for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
      RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
      I1 = iI1*7;
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termCode.cc>
#if 0
      if (iI1==0) {
        plabmx = pxm;
        plabmy = pym;
        plabmz = pzm;
      } else {
        plabmx = 0.0;
        plabmy = 0.0;
        plabmz = 0.0;
      }
#endif
      parts->vectorPushExtend(plabmx);
      parts->vectorPushExtend(plabmy);
      parts->vectorPushExtend(plabmz);
#if 0
      if (iI1<3) {
        core::write_bf_stream(BF("%s:%d  I1cur[%d] I1[%d] -> Pos         %f %f %f\n") % __FILE__ % __LINE__ % I1cur % I1 % ea1._Position.getX() % ea1._Position.getY() % ea1._Position.getZ() );
        core::write_bf_stream(BF("%s:%d  I1cur[%d] I1[%d] -> pxm,pym,pzm %f %f %f\n") % __FILE__ % __LINE__ % I1cur % I1 % pxm % pym % pzm );
        core::write_bf_stream(BF("%s:%d               point -> %f %f %f\n") % __FILE__ % __LINE__ % plabmx % plabmy %plabmz );
      }
#endif
    }
    I1start = I1end;
  }
  return parts;
}


CL_DEFMETHOD core::ComplexVector_sp EnergyRigidBodyNonbond_O::write_nonbond_atoms_to_complex_vector(core::ComplexVector_sp parts)
{
  size_t istart = 0;
#undef	NONBOND_POSITION_RB_SET_PARAMETER
#define	NONBOND_POSITION_RB_SET_PARAMETER(x)	{}
#undef	NONBOND_POSITION_RB_SET_POSITION
#define	NONBOND_POSITION_RB_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_POSITION_RB_SET_POINT
#define	NONBOND_POSITION_RB_SET_POINT(x,ii,of)	{x=ii._Position.of;}
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termDeclares.cc>
#pragma clang diagnostic pop
  double am, bm, cm, dm, xm, ym, zm;
  double pxm, pym, pzm;
  int	I1;
  size_t I1start = 0;
  for ( size_t iI1 = 0; iI1<this->_RigidBodyEndAtom->length(); ++iI1 ) {
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
    for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
      RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
      parts->vectorPushExtend(ea1._Object);
    }
    I1start = I1end;
  }
  return parts;
}


void	EnergyRigidBodyNonbond_O::setupHessianPreconditioner(
                                                    NVector_sp nvPosition,
                                                    AbstractLargeSquareMatrix_sp m )
{
  SIMPLE_ERROR(BF("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!!"));
}


inline double periodic_boundary_adjust(const double& delta, const double& rsize, const double& size)
{
  double result = delta;
  result -= size*std::nearbyint(delta*rsize);
  return result;
}

double	EnergyRigidBodyNonbond_O::evaluateAllComponent( ScoringFunction_sp score,
                                                        NVector_sp 	pos,
                                                        bool 		calcForce,
                                                        gc::Nilable<NVector_sp> 	force,
                                                        bool		calcDiagonalHessian,
                                                        bool		calcOffDiagonalHessian,
                                                        gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                        gc::Nilable<NVector_sp>	hdvec, 
                                                        gc::Nilable<NVector_sp> 	dvec )
{
  if (this->_CrossTerms.size() == 0 ) this->initializeCrossTerms(false);
  double vdwScale = this->getVdwScale();
  double electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
  this->_EnergyElectrostatic = 0.0;
  this->_EnergyVdw = 0.0;
  this->_TotalEnergy = 0.0;
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
  RigidBodyEnergyFunction_sp rigidBodyEnergyFunction = gc::As<RigidBodyEnergyFunction_sp>(score);
  BoundingBox_sp boundingBox = rigidBodyEnergyFunction->boundingBox();
  if (boundingBox.unboundp()) {
    SIMPLE_ERROR(BF("The rigid-body-energy-function bounding-box is unbound - it must be defined"));
  }
  Vector3 rwidths = boundingBox->get_bounding_box_rwidths();
  double x_rsize = rwidths.getX();
  double y_rsize = rwidths.getY();
  double z_rsize = rwidths.getZ();
  Vector3 widths = boundingBox->get_bounding_box_widths();
  double x_size = widths.getX();
  double y_size = widths.getY();
  double z_size = widths.getZ();
  double half_x_size = x_size/2.0;
  double half_y_size = y_size/2.0;
  double half_z_size = z_size/2.0;
  
#define NONBONDRB_CALC_FORCE
#define NONBONDRB_CALC_DIAGONAL_HESSIAN
#define NONBONDRB_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBONDRB_SET_PARAMETER
#define	NONBONDRB_SET_PARAMETER(x)	{}
#undef	NONBONDRB_SET_POSITION
#define	NONBONDRB_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBONDRB_SET_POINT
#define	NONBONDRB_SET_POINT(x,ii,of)	{x=ii._Position.of;}
#undef	NONBONDRB_EEEL_ENERGY_ACCUMULATE
#define	NONBONDRB_EEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);} // not used
#undef	NONBONDRB_EVDW_ENERGY_ACCUMULATE
#define	NONBONDRB_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);} // not used
#undef	NONBONDRB_ENERGY_ACCUMULATE
#define	NONBONDRB_ENERGY_ACCUMULATE(e) {this->_TotalEnergy += e;};
#undef	NONBONDRB_FORCE_ACCUMULATE
#undef	NONBONDRB_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBONDRB_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBONDRB_FORCE_ACCUMULATE 		ForceAcc
#define	NONBONDRB_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBONDRB_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
#define PBX(_delta_) periodic_boundary_adjust(_delta_,x_rsize,x_size)
#define PBY(_delta_) periodic_boundary_adjust(_delta_,y_rsize,y_size)
#define PBZ(_delta_) periodic_boundary_adjust(_delta_,z_rsize,z_size)
#define BAIL_OUT_IF_CUTOFF(deltax,deltay,deltaz) \
  if (fabs(deltax)>half_x_size) continue; \
  if (fabs(deltay)>half_y_size) continue; \
  if (fabs(deltaz)>half_z_size) continue;
  
  if ( !this->isEnabled() ) return 0.0;
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
  LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_termDeclares.cc>
#pragma clang diagnostic pop
  double dA,dC,dQ1Q2;
  double am, bm, cm, dm, xm, ym, zm;
  double pxm, pym, pzm;
  double an, bn, cn, dn, xn, yn, zn;
  double pxn, pyn, pzn;
  int	I1, I2;
  size_t interactions = 0;
  size_t I1start = 0;
  for ( size_t iHelix1 = 0; iHelix1<(this->_RigidBodyEndAtom->length()-1); ++iHelix1 ) {
    size_t I1end = (*this->_RigidBodyEndAtom)[iHelix1];
    for (size_t iHelix2 = iHelix1+1; iHelix2<this->_RigidBodyEndAtom->length(); ++iHelix2 ) {
      size_t I2start = (*this->_RigidBodyEndAtom)[iHelix2-1];
      size_t I2end = (*this->_RigidBodyEndAtom)[iHelix2];
      for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
        RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
        double charge1 = ea1._Charge;
        double electrostatic_scaled_charge1 = charge1*electrostaticScale;
        for ( size_t I2cur = I2start; I2cur<I2end; ++I2cur ) {
          ++interactions;
          RigidBodyAtomInfo& ea2 = this->_AtomInfoTable[I2cur];
          RigidBodyNonbondCrossTerm& crossTerm = this->crossTerm(ea1._TypeIndex,ea2._TypeIndex);
          dA = crossTerm.dA;
          dC = crossTerm.dC;
#if 0          
          double rStar = ea1._Radius+ea2._Radius;
          double epsilonij = sqrt(ea1._Epsilon*ea2._Epsilon);
          double rStar2 = rStar*rStar;
          double rStar6 = rStar2*rStar2*rStar2;
          double rStar12 = rStar6*rStar6;
          dA = epsilonij*rStar12*vdwScale;
          dC = 2.0*epsilonij*rStar6*vdwScale;
#endif
          double charge2 = ea2._Charge;
          dQ1Q2 = electrostatic_scaled_charge1*charge2;
          I1 = iHelix1*7;
          I2 = iHelix2*7;
#include <cando/chem/energy_functions/_NONBONDRBPB_termCode.cc>
#if 0
          if (Energy > 1.0) {
            core::write_bf_stream(BF("Energy iHelix1 iHelix2 I1cur I2cur\n"));
            core::write_bf_stream(BF("Energy= %f %d %d %d %d\n") % Energy % iHelix1 % iHelix2 %I1cur % I2cur);
            core::write_bf_stream(BF("am,bm,cm,dm,xm,ym,zm,pxm,pym,pzm -> %f,%f,%f,%f,%f,%f,%f %f,%f,%f\n") % am % bm % cm % dm % xm % ym % zm % pxm % pym % pzm);
            core::write_bf_stream(BF("an,bn,cn,dn,xn,yn,zn,pxn,pyn,pzn -> %f,%f,%f,%f,%f,%f,%f %f,%f,%f\n") % an % bn % cn % dn % xn % yn % zn % pxn % pyn % pzn);
            core::write_bf_stream(BF("DeltaX, DeltaY, DeltaZ  -> %f, %f, %f\n") % DeltaX % DeltaY % DeltaZ );
            core::write_bf_stream(BF("dA, dC -> %f, %f\n") % dA % dC );
          }
          
          if (std::isnan(Energy)) {
            SIMPLE_ERROR(BF("Energy is nan"));
          }
#endif
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
          nbi->_calcForce = calcForce;
          nbi->_calcDiagonalHessian = calcDiagonalHessian;
          nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/chem/energy_functions/_NONBONDRBPB_debugEvalSet.cc>
#endif //]
#ifdef DEBUG_NONBOND_TERM
          if ( this->_DebugEnergy ) {
            std::string key;
            key = _rep_(ea1._Object)+"-"+_rep_(ea2._Object);
            LOG_ENERGY(BF( "MEISTER nonbond %s args cando\n")% key );
            LOG_ENERGY(BF( "MEISTER nonbond %s iHelix1->%3d iHelix2->%3d  I1Cur->%3d I2Cur->%3d\n")% key % iHelix1 % iHelix2 % I1cur % I2cur );
            LOG_ENERGY(BF( "MEISTER nonbond %s dA %5.3lf\n")% key % dA );
            LOG_ENERGY(BF( "MEISTER nonbond %s dC %5.3lf\n")% key % dC );
            LOG_ENERGY(BF( "MEISTER nonbond %s dQ1Q2 %5.3lf\n")% key % dQ1Q2 );
#if 0
            LOG_ENERGY(BF( "MEISTER nonbond %s x1 %5.3lf %d\n")% key % x1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s y1 %5.3lf %d\n")% key % y1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s z1 %5.3lf %d\n")% key % z1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s x2 %5.3lf %d\n")% key % x2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s y2 %5.3lf %d\n")% key % y2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s z2 %5.3lf %d\n")% key % z2 % (I2/3+1) );
#endif
            LOG_ENERGY(BF( "MEISTER nonbond %s results\n")% key );
            LOG_ENERGY(BF( "MEISTER nonbond %s EnergyVdw %lf\n")% key % EnergyVdw);
            LOG_ENERGY(BF( "MEISTER nonbond %s EnergyElectrostatic %lf\n")% key % EnergyElectrostatic);
            LOG_ENERGY(BF( "MEISTER nonbond %s Enonbond(EnergyVdw+EnergyElectrostatic) %lf\n")% key % (EnergyVdw+EnergyElectrostatic) );
            if ( calcForce ) {
              LOG_ENERGY(BF( "MEISTER nonbond %s fam %lf %d\n")% key % fam % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fbm %lf %d\n")% key % fbm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fcm %lf %d\n")% key % fcm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fdm %lf %d\n")% key % fdm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fxm %lf %d\n")% key % fxm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fym %lf %d\n")% key % fym % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fzm %lf %d\n")% key % fzm % (I1cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fan %lf %d\n")% key % fan % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fbn %lf %d\n")% key % fbn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fcn %lf %d\n")% key % fcn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fdn %lf %d\n")% key % fdn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fxn %lf %d\n")% key % fxn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fyn %lf %d\n")% key % fyn % (I2cur/3+1) );
              LOG_ENERGY(BF( "MEISTER nonbond %s fzn %lf %d\n")% key % fzn % (I2cur/3+1) );
            }
            LOG_ENERGY(BF( "MEISTER nonbond %s stop\n")% key );
          }
#endif
        }
      }
    }
    gctools::handle_all_queued_interrupts();
    I1start = I1end;
  }
#ifdef DEBUG_NONBOND_TERM
  if ( this->_DebugEnergy ) {
    LOG_ENERGY(BF( "MEISTER nonbond interactions -> %d\n") % interactions);
  }
#endif
  
#if 0
  printf("%s:%d NonbondRigidBody    TotalEnergy -> %lf\n", __FILE__, __LINE__, this->_TotalEnergy );
#endif
  if (std::isnan(this->_TotalEnergy)) {
    SIMPLE_ERROR(BF("Returning nan"));
  }
#if 0  
  core::write_bf_stream(BF("Total-energy %e\n") % this->_TotalEnergy );
#endif
  return this->_TotalEnergy;
}
    
void	EnergyRigidBodyNonbond_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                                NVector_sp 	pos)
{_OF();
  IMPLEMENT_ME();
  // Use the code below
#if 0
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//------------------


#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{x=nbi->term.x;}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE(i,o,v) {}
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}

  if ( this->isEnabled() ) {
    _BLOCK_TRACE("NonbondEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_termDeclares.cc>
#pragma clang diagnostic pop
    double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
    int	I1, I2,i;
    gctools::Vec0<EnergyRigidBodyNonbond>::iterator nbi;
    for ( i=0,nbi=this->_Terms.begin();
          nbi!=this->_Terms.end(); nbi++,i++ ) {
#include <cando/chem/energy_functions/_NONBONDRBPB_termCode.cc>
      int index = i;
#include <cando/chem/energy_functions/_NONBONDRBPB_debugFiniteDifference.cc>

    }
  }
#endif
}




void EnergyRigidBodyNonbond_O::initialize()
{
  this->Base::initialize();
  this->setErrorThreshold(1.0);
  this->setDielectricConstant(80.0);
  this->setVdwScale(1.0);
  this->setElectrostaticScale(1.0);
}






};
