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
HE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
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
#include <cando/geom/ovector3.h>
#include <cando/chem/loop.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

#if 0
#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
#define LOG_ENERGY BFORMAT_T
#endif

namespace chem
{

core::List_sp RigidBodyNonbondCrossTerm::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,da),mk_double_float(this->dA)), 
                                  core::Cons_O::create(INTERN_(kw,dc),mk_double_float(this->dC)) );
}

void RigidBodyNonbondCrossTerm::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of RigidBodyNonbondCrossTerm");
}


core::List_sp RigidBodyAtomInfo::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,object),this->_Object), 
                                  core::Cons_O::create(INTERN_(kw,typeIndex),core::make_fixnum(this->_TypeIndex)),
                                  core::Cons_O::create(INTERN_(kw,radius), mk_double_float(this->_Radius)),
                                  core::Cons_O::create(INTERN_(kw,epsilon), mk_double_float(this->_Epsilon)),
                                  core::Cons_O::create(INTERN_(kw,charge), mk_double_float(this->_Charge)),
                                  core::Cons_O::create(INTERN_(kw,position), translate::to_object<Vector3>::convert(this->_Position) ) );
}

void RigidBodyAtomInfo::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of RigidBodyAtomInfo");
}

struct TypeParameters {
  num_real _Radius;
  num_real _Epsilon;
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
    core::clasp_write_string(fmt::format("There are {} types\n" , this->_NumberOfTypes));
  }
  this->_CrossTerms.resize(next_type_index*next_type_index);
  num_real vdwScale = this->getVdwScale();
  for ( size_t xi=0; xi<next_type_index; xi++ ) {
    for (size_t yi=xi; yi<next_type_index; yi++ ) {
      TypeParameters& ea1 = typeParameters[xi];
      TypeParameters& ea2 = typeParameters[yi];
      num_real rStar = ea1._Radius+ea2._Radius;
      num_real epsilonij = sqrt(ea1._Epsilon*ea2._Epsilon);
      num_real rStar2 = rStar*rStar;
      num_real rStar6 = rStar2*rStar2*rStar2;
      num_real rStar12 = rStar6*rStar6;
      num_real dA = epsilonij*rStar12*vdwScale;
      num_real dC = 2.0*epsilonij*rStar6*vdwScale;
      RigidBodyNonbondCrossTerm crossTerm(dA,dC);
      if (verbose) {
        core::clasp_write_string(fmt::format("Cross term for types {} - {}\n" , xi , yi ));
        core::clasp_write_string(fmt::format("    type: {}   radius-> {}   epsilon-> {}\n" , xi , ea1._Radius , ea1._Epsilon));
        core::clasp_write_string(fmt::format("    type: {}   radius-> {}   epsilon-> {}\n" , yi , ea2._Radius , ea2._Epsilon));
        core::clasp_write_string(fmt::format("    dA -> {}      dC -> {}\n" , dA , dC ));
      }
      this->crossTerm(xi,yi) = crossTerm;
      this->crossTerm(yi,xi) = crossTerm;
    }
  }
}
          
EnergyRigidBodyNonbond_sp EnergyRigidBodyNonbond_O::make(core::Array_sp end_atoms) {
  if (end_atoms->length()<1) {
    SIMPLE_ERROR("You must provide a vector of end atom indexes with at least one end atom");
  }
  core::SimpleVector_byte32_t_sp sv = core::SimpleVector_byte32_t_O::make(end_atoms->length());
  // Copy the atom indexes
  size_t istart = 0;
  for ( size_t i(0), iEnd(end_atoms->length()); i<iEnd; ++i ) {
    Fixnum f = core::clasp_to_fixnum(end_atoms->rowMajorAref(i));
    if (istart>=f) {
      SIMPLE_ERROR("The list of atom indexes must be strictly increasing");
    }
    (*sv)[i] = f;
  }
  auto nb = gctools::GC<EnergyRigidBodyNonbond_O>::allocate(sv);
  return nb;
}

 void EnergyRigidBodyNonbond_O::energyRigidBodyNonbondSetTerm(gc::Fixnum index, core::T_sp atom, double radius, double epsilon, double charge, const Vector3& position) {
  if (index < 0 || index >= this->_AtomInfoTable.size()) {
    SIMPLE_ERROR("Index out of range {} - max is {}" , index , this->_AtomInfoTable.size());
  }
  if (this->_AtomInfoTable[index]._Object.boundp()) {
    SIMPLE_ERROR("The rigid body nonbond term at {} has already been set to object {}" , index , _rep_(this->_AtomInfoTable[index]._Object));
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
  SIMPLE_ERROR("index {} is out of bounds" , index);
}


void	EnergyRigidBodyNonbond_O::dumpTerms(core::HashTable_sp atomTypes)
{
  size_t istart = 0;
  for ( size_t rb = 0; rb<this->_RigidBodyEndAtom->length(); ++rb) {
    for ( size_t i = istart; i<(*this->_RigidBodyEndAtom)[rb]; ++i ) {
      RigidBodyAtomInfo& ai = this->_AtomInfoTable[i];
      core::clasp_write_string(fmt::format("I1 = {:3}  {}  {:f} {:f} {:f} {}\n" , (rb*7) , _rep_(ai._Object) , ai._Radius , ai._Epsilon , ai._Charge , ai._Position.asString()));
    }
    istart = (*this->_RigidBodyEndAtom)[rb];
  }
}

CL_DEFMETHOD core::List_sp EnergyRigidBodyNonbond_O::parts_as_list(NVector_sp pos)
{
  ql::list result;
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
  num_real am, bm, cm, dm, xm, ym, zm;
  num_real pxm, pym, pzm;
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

size_t EnergyRigidBodyNonbond_O::partsCoordinates(NVector_sp pos, size_t idx, core::SimpleVector_float_sp coords )
{
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
  num_real am, bm, cm, dm, xm, ym, zm;
  num_real pxm, pym, pzm;
  int	I1;
  size_t I1start = 0;
  for ( size_t iI1 = 0; iI1<this->_RigidBodyEndAtom->length(); ++iI1 ) {
    ql::list helix;
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
    for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
      RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
      I1 = iI1*7;
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termCode.cc>
      (*coords)[idx++] = plabmx;
      (*coords)[idx++] = plabmy;
      (*coords)[idx++] = plabmz;
    }
    I1start = I1end;
  }
  return idx;
}


struct coordinate_lookup {
  coordinate_lookup(NVector_sp coords) : _Coordinates(coords) {};
  NVector_sp _Coordinates;
  size_t _StartIndex;
  void setStartIndex(size_t startIndex) { this->_StartIndex = startIndex; };
  num_real getX() const { return (*this->_Coordinates)[this->_StartIndex+0];}
  num_real getY() const { return (*this->_Coordinates)[this->_StartIndex+1];}
  num_real getZ() const { return (*this->_Coordinates)[this->_StartIndex+2];}
};



CL_DOCSTRING(R"dx(This function is to test the rigid body transformation code.
Arguments:
rigid-body-pos : An nvector of N quaternion/origin coordinates (a,b,c,d,x,y,z).
end-indexesx3 : A vector of integer indexes that index into 'coordinates' and indicate where
each object (defined by a contiguous set of points) in coordinates ends. The 'end-indexesx3'
values are indexes directly into 'coordinates', so they are x3.
coordinates : An nvector of coordinates.
output : A complex-vector-float where the transformed points are written.
)dx")
CL_DEFMETHOD core::ComplexVector_float_sp EnergyRigidBodyNonbond_O::write_rigid_body_coordinates_to_complex_vector_float(NVector_sp rigid_body_pos, core::Array_sp end_indexesx3, NVector_sp coordinates, core::ComplexVector_float_sp output)
{
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
  num_real am, bm, cm, dm, xm, ym, zm;
  num_real pxm, pym, pzm;
  int	I1;
  size_t I1start = 0;
  coordinate_lookup ea1(coordinates); 
  for ( size_t iI1 = 0; iI1<end_indexesx3->length(); ++iI1 ) {
    size_t I1end = core::clasp_to_size(end_indexesx3->rowMajorAref(iI1));
    for ( size_t I1cur = I1start; I1cur<I1end; I1cur+=3 ) {
      ea1.setStartIndex(I1cur);
      I1 = iI1*7;
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termCode.cc>
      output->vectorPushExtend((float)plabmx);
      output->vectorPushExtend((float)plabmy);
      output->vectorPushExtend((float)plabmz);
#if 0
      if (coordIndex==0) {
        core::clasp_write_string(fmt::format("fill_pointer -> {}\n" , output->fillPointer() ));
        core::clasp_write_string(fmt::format("&plabmx -> {}\n" , (void*)&plabmx));
        core::clasp_write_string(fmt::format("&output[0] -> {}\n" , (void*)&(*output)[0]));
        core::clasp_write_string(fmt::format("Point[0] -> {}, {}, {}\n" , plabmx , plabmy , plabmz));
        core::clasp_write_string(fmt::format("output[0] -> {}\n" , (*output)[0]));
        core::clasp_write_string(fmt::format("hex dump output[0] -> {}X\n" , *(unsigned int*)&(*output)[0]));
        (*output)[0] = plabmx;
        core::clasp_write_string(fmt::format("after output[0] -> {}\n" , (*output)[0]));
        core::clasp_write_string(fmt::format("hex dump output[0] -> {}X\n" , *(unsigned int*)&(*output)[0]));
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
  num_real am, bm, cm, dm, xm, ym, zm;
  num_real pxm, pym, pzm;
  int	I1;
  size_t I1start = 0;
  for ( size_t iI1 = 0; iI1<this->_RigidBodyEndAtom->length(); ++iI1 ) {
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
#if 0
    if (iI1<3) {
      core::clasp_write_string(fmt::format("{}:{} I1start -> {}  I1end -> {}\n" , __FILE__ , __LINE__ , I1start , I1end));
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
        core::clasp_write_string(fmt::format("{}:{}  I1cur[{}] I1[{}] -> Pos         {} {} {}\n" , __FILE__ , __LINE__ , I1cur , I1 , ea1._Position.getX() , ea1._Position.getY() , ea1._Position.getZ() ));
        core::clasp_write_string(fmt::format("{}:{}  I1cur[{}] I1[{}] -> pxm,pym,pzm {} {} {}\n" , __FILE__ , __LINE__ , I1cur , I1 , pxm , pym , pzm ));
        core::clasp_write_string(fmt::format("{}:{}               point -> {} {} {}\n" , __FILE__ , __LINE__ , plabmx , plabmy {}labmz ));
      }
#endif
    }
    I1start = I1end;
  }
  return parts;
}


CL_DEFMETHOD core::ComplexVector_sp EnergyRigidBodyNonbond_O::write_nonbond_atoms_to_complex_vector(core::ComplexVector_sp parts)
{
#undef	NONBOND_POSITION_RB_SET_PARAMETER
#define	NONBOND_POSITION_RB_SET_PARAMETER(x)	{}
#undef	NONBOND_POSITION_RB_SET_POSITION
#define	NONBOND_POSITION_RB_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_POSITION_RB_SET_POINT
#define	NONBOND_POSITION_RB_SET_POINT(x,ii,of)	{x=ii._Position.of;}
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_POSITIONS_termDeclares.cc>
  num_real am, bm, cm, dm, xm, ym, zm;
  num_real pxm, pym, pzm;
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
#pragma clang diagnostic pop
  return parts;
}


void	EnergyRigidBodyNonbond_O::setupHessianPreconditioner(
                                                    NVector_sp nvPosition,
                                                    AbstractLargeSquareMatrix_sp m,
                                                    core::T_sp activeAtomMask )
{
  SIMPLE_ERROR("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!!");
}


inline num_real periodic_boundary_adjust(const num_real& delta, const num_real& rsize, const num_real& size)
{
  num_real result = delta;
  result -= size*std::nearbyint(delta*rsize);
  return result;
}

double EnergyRigidBodyNonbond_O::evaluateAllComponent( ScoringFunction_sp score,
                                                       NVector_sp 	pos,
                                                       core::T_sp energyScale,
                                                       core::T_sp energyComponents,
                                                       bool 		calcForce,
                                                       gc::Nilable<NVector_sp> 	force,
                                                       bool		calcDiagonalHessian,
                                                       bool		calcOffDiagonalHessian,
                                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                       gc::Nilable<NVector_sp>	hdvec, 
                                                       gc::Nilable<NVector_sp> 	dvec,
                                                       core::T_sp activeAtomMask,
                                                       core::T_sp debugInteractions )
{
  double dielectricConstant;
  double dQ1Q2Scale;
  double cutoff;
  energyFunctionNonbondParameters(score,energyScale,dielectricConstant,dQ1Q2Scale,cutoff);
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
//  SIMPLE_WARN("FIXactiveAtomMask How do I deal with activeAtomMask");
  this->_Evaluations++;
  if (this->_CrossTerms.size() == 0 ) this->initializeCrossTerms(false);
  num_real electrostaticScale = this->getElectrostaticScale()*dQ1Q2Scale/dielectricConstant;
  double totalEnergy = 0.0;
  bool	hasForce = force.notnilp();
  RigidBodyEnergyFunction_sp rigidBodyEnergyFunction = gc::As<RigidBodyEnergyFunction_sp>(score);
  BoundingBox_sp boundingBox = rigidBodyEnergyFunction->boundingBox();
  if (boundingBox.unboundp()) {
    SIMPLE_ERROR("The rigid-body-energy-function bounding-box is unbound - it must be defined");
  }
  Vector3 rwidths = boundingBox->get_bounding_box_rwidths();
  num_real x_rsize = rwidths.getX();
  num_real y_rsize = rwidths.getY();
  num_real z_rsize = rwidths.getZ();
  Vector3 widths = boundingBox->get_bounding_box_widths();
  num_real x_size = widths.getX();
  num_real y_size = widths.getY();
  num_real z_size = widths.getZ();
  num_real half_x_size = x_size/2.0;
  num_real half_y_size = y_size/2.0;
  num_real half_z_size = z_size/2.0;
  
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
#define	NONBONDRB_EEEL_ENERGY_ACCUMULATE(e) {energyElectrostatic +=(e);} // not used
#undef	NONBONDRB_EVDW_ENERGY_ACCUMULATE
#define	NONBONDRB_EVDW_ENERGY_ACCUMULATE(e) {energyVdw+=(e);} // not used
#undef	NONBONDRB_ENERGY_ACCUMULATE
#define	NONBONDRB_ENERGY_ACCUMULATE(e) {totalEnergy += e;};
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
  
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
  LOG("Nonbond component is enabled" );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_termDeclares.cc>
#pragma clang diagnostic pop
  num_real dA,dC,dQ1Q2;
  num_real am, bm, cm, dm, xm, ym, zm;
  num_real pxm, pym, pzm;
  num_real an, bn, cn, dn, xn, yn, zn;
  num_real pxn, pyn, pzn;
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
        num_real charge1 = ea1._Charge;
        num_real electrostatic_scaled_charge1 = charge1*electrostaticScale;
        for ( size_t I2cur = I2start; I2cur<I2end; ++I2cur ) {
          ++interactions;
          RigidBodyAtomInfo& ea2 = this->_AtomInfoTable[I2cur];
          RigidBodyNonbondCrossTerm& crossTerm = this->crossTerm(ea1._TypeIndex,ea2._TypeIndex);
          dA = crossTerm.dA;
          dC = crossTerm.dC;
#if 0          
          num_real rStar = ea1._Radius+ea2._Radius;
          num_real epsilonij = sqrt(ea1._Epsilon*ea2._Epsilon);
          num_real rStar2 = rStar*rStar;
          num_real rStar6 = rStar2*rStar2*rStar2;
          num_real rStar12 = rStar6*rStar6;
          dA = epsilonij*rStar12*vdwScale;
          dC = 2.0*epsilonij*rStar6*vdwScale;
#endif
          num_real charge2 = ea2._Charge;
          dQ1Q2 = electrostatic_scaled_charge1*charge2;
          I1 = iHelix1*7;
          I2 = iHelix2*7;
#include <cando/chem/energy_functions/_NONBONDRBPB_termCode.cc>
#if 0
          if (Energy > 1.0) {
            core::clasp_write_string("Energy iHelix1 iHelix2 I1cur I2cur\n");
            core::clasp_write_string(fmt::format("Energy= {} {} {} {} {}\n" , Energy , iHelix1 , iHelix2, I1cur, I2cur));
            core::clasp_write_string(fmt::format("am,bm,cm,dm,xm,ym,zm,pxm,pym,pzm -> {},{},{},{},{},{},{} {},{},{}\n" , am , bm , cm , dm , xm , ym , zm , pxm , pym , pzm));
            core::clasp_write_string(fmt::format("an,bn,cn,dn,xn,yn,zn,pxn,pyn,pzn -> {},{},{},{},{},{},{} {},{},{}\n" , an , bn , cn , dn , xn , yn , zn , pxn , pyn , pzn));
            core::clasp_write_string(fmt::format("DeltaX, DeltaY, DeltaZ  -> {}, {}, {}\n" , DeltaX , DeltaY , DeltaZ ));
            core::clasp_write_string(fmt::format("dA, dC -> {}, {}\n" , dA , dC ));
          }
          
          if (std::isnan(Energy)) {
            SIMPLE_ERROR("Energy is nan");
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
        }
      }
    }
    gctools::handle_all_queued_interrupts();
    I1start = I1end;
  }
  maybeSetEnergy( energyComponents, EnergyRigidBodyNonbond_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}
    
void	EnergyRigidBodyNonbond_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                                NVector_sp 	pos)
{
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

  {
    
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_NONBONDRBPB_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
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
  this->setDielectricConstant(80.0);
  this->setVdwScale(1.0);
  this->setElectrostaticScale(1.0);
}


void EnergyRigidBodyNonbond_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,dielectricConstant) ,this->_DielectricConstant );
  node->field( INTERN_(kw,scaleVdw), this->_ScaleVdw );
  node->field( INTERN_(kw,scaleElectrostatic), this->_ScaleElectrostatic );
  node->field( INTERN_(kw,energyVdw), this->_EnergyVdw );
  node->field( INTERN_(kw,energyElectrostatic), this->_EnergyElectrostatic);
  node->field( INTERN_(kw,rigidBodyEndAtom), this->_RigidBodyEndAtom );
  node->field( INTERN_(kw,atomInfoTable), this->_AtomInfoTable );
  node->field( INTERN_(kw,numberOfTypes), this->_NumberOfTypes );
  node->field( INTERN_(kw,crossTerms), this->_CrossTerms );
  node->Base::fields(node);
}

};
