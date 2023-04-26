//#define DEBUG_LEVEL_FULL
/*
    File: energyAtomTable.cc
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

//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/foundation.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/nativeVector.h>
#include <clasp/core/ql.h>
#include <clasp/core/bformat.h>
#include <clasp/core/sort.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/array.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/atom.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/residue.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/forceField.h>
#include <cando/units/unitsPackage.h>
#include <cando/units/quantity.h>
#include <cando/chem/loop.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/ffNonbondDb.h>
#include <clasp/core/translators.h>
#include <clasp/core/wrappers.h>




namespace chem
{
EnergyAtom::EnergyAtom() : _Flag(0) {};
EnergyAtom::EnergyAtom(Atom_sp atom,uint coordinateIndex) : Base(atom,coordinateIndex), _Flag(0)  {} ;
EnergyAtom::EnergyAtom(core::T_sp forceField, Atom_sp atom, uint coordinateIndex, core::HashTable_sp atomTypes ) : Base(atom,coordinateIndex), _Flag(0)
{
  this->defineForAtom(forceField,atom,coordinateIndex,atomTypes);
}


void EnergyAtom::defineForAtom(core::T_sp forceField, Atom_sp a1, uint coordinateIndex, core::HashTable_sp atomTypes )
{
  this->setupBase(a1,coordinateIndex);
  if (a1->getType(atomTypes).nilp()) {
    SIMPLE_ERROR(("The atom type of %s is NIL!") , _rep_(a1));
  }
  core::T_sp typeIndex = core::eval::funcall(_sym_find_atom_type_position,forceField,a1->getType(atomTypes));
  if (!typeIndex.fixnump()) {
    TYPE_ERROR(typeIndex,cl::_sym_fixnum);
  }
  this->_TypeIndex = typeIndex.unsafe_fixnum();
  this->_Charge = a1->getCharge();
  this->_AtomicNumber = a1->getAtomicNumber();
}

#ifdef XML_ARCHIVE
void	EnergyAtom::archive(core::ArchiveP node)
{
  this->Base::archive(node);
}
#endif


#if 0 //[
adapt::QDomNode_sp	EnergyAtom::asXmlRelativeToContainer(chem::Matter_sp parent)
{
  adapt::QDomNode_sp	node;
  string	path;
  node = chem::new_adapt::QDomNode_sp("EnergyAtom");
  path = this->_Atom->getIdPath(parent);
  node->addAttributeString("storagePath",path);
  node->addAttributeString("atomName",this->_Atom->getName());
  node->addAttributeInt("coordinateIndex",this->_CoordinateIndex);
  node->addAttributeDoubleScientific("x",this->_Atom->getPosition().getX());
  node->addAttributeDoubleScientific("y",this->_Atom->getPosition().getY());
  node->addAttributeDoubleScientific("z",this->_Atom->getPosition().getZ());
  node->addAttributeDoubleScientific("fx",this->_Atom->getForce().getX());
  node->addAttributeDoubleScientific("fy",this->_Atom->getForce().getY());
  node->addAttributeDoubleScientific("fz",this->_Atom->getForce().getZ());
  return node;
}


void	EnergyAtom::parseFromXmlRelativeToContainer(adapt::QDomNode_sp xml,
                                                    chem::Matter_sp parent)
{
  string	path;
  double	x,y,z,fx,fy,fz;
  Vector3	v;
  LOG("parsing xml with localName=%s" , xml->getLocalName().c_str() );
  path = xml->getAttributeString("storagePath");
  this->_Atom = (parent->getContentForIdPath(path)).as<chem::Atom_O>();
  ASSERTNOTNULL(this->_Atom);
  this->_CoordinateIndex = xml->getAttributeInt("coordinateIndex");
  x = xml->getAttributeDouble("x");
  y = xml->getAttributeDouble("y");
  z = xml->getAttributeDouble("z");
  v.set(x,y,z);
  this->_Atom->setPosition(v);
  fx = xml->getAttributeDouble("fx");
  fy = xml->getAttributeDouble("fy");
  fz = xml->getAttributeDouble("fz");
  v.set(fx,fy,fz);
  this->_Atom->setForce(v);
}
#endif //]

string		EnergyAtom::getResidueAndName(Residue_sp res)
{
  stringstream	ss;
  ss.str("");
  ss << ":"<<res->getId()<<"@"<<this->atom()->getName();
  return ss.str();
}

CL_LISPIFY_NAME(make-atom-table);
CL_DEF_CLASS_METHOD
AtomTable_sp AtomTable_O::make()
{
  auto eat = gctools::GC<AtomTable_O>::allocate_with_default_constructor();
  return eat;
}

void AtomTable_O::initialize()
{
  this->_AtomTableIndices = core::HashTableEq_O::create_default();
  core::ComplexVector_int32_t_sp residue_pointers = core::ComplexVector_int32_t_O::make_vector(32,0,core::make_fixnum(0),nil<core::T_O>(),false,core::make_fixnum(0));
  core::ComplexVector_T_sp residue_names = core::ComplexVector_T_O::make(32,nil<core::T_O>(),core::make_fixnum(0));
  core::ComplexVector_int32_t_sp atoms_per_molecule = core::ComplexVector_int32_t_O::make_vector(32,0,core::make_fixnum(0),nil<core::T_O>(),false,core::make_fixnum(0));
//  atoms_per_molecule->vectorPushExtend(0);
  this->_ResiduePointers = residue_pointers;
  this->_ResidueNames = residue_names;
  this->_AtomsPerMolecule = atoms_per_molecule;
  this->_Residues = core::core__make_vector(_lisp->_true(), 16, true, core::make_fixnum(0));
  this->_Molecules = core::core__make_vector(_lisp->_true(), 16, true, core::make_fixnum(0));
}

core::List_sp EnergyAtom::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,atom_name),this->_AtomName),
                                  core::Cons_O::create(INTERN_(kw,charge),core::clasp_make_double_float(this->_Charge)),
                                  core::Cons_O::create(INTERN_(kw,mass),core::clasp_make_double_float(this->_Mass)),
                                  core::Cons_O::create(INTERN_(kw,type_index), core::make_fixnum(this->_TypeIndex)),
                                  core::Cons_O::create(INTERN_(kw,atomic_number), core::make_fixnum(this->_AtomicNumber)));
}



void AtomTable_O::fields(core::Record_sp node) {
  node->field( INTERN_(kw,atoms), this->_Atoms );
  node->field( INTERN_(kw,atom_table_indices), this->_AtomTableIndices);
  node->field( INTERN_(kw,residue_pointers), this->_ResiduePointers);
  node->field( INTERN_(kw,residue_names), this->_ResidueNames);
  node->field( INTERN_(kw,atoms_per_molecule),this->_AtomsPerMolecule);
  node->field( INTERN_(kw,residues),this->_Residues);
  node->field( INTERN_(kw,molecules),this->_Molecules);
  this->Base::fields(node);
}

CL_DEFMETHOD size_t AtomTable_O::getCoordinateIndex(Atom_sp a)
{
  EnergyAtom* ea = this->getEnergyAtomPointer(a);
  return ea->coordinateIndexTimes3();
}

CL_DOCSTRING("Read coordinates from atoms into an NVector of coordinates");
CL_DEFMETHOD void AtomTable_O::readAtomCoordinates(NVector_sp coords) const
{
  if (coords->length() != this->getNVectorSize()) {
    SIMPLE_ERROR("The passed NVector has length %d but the atom-table needs length %d", coords->length(), this->getNVectorSize() );
  }
  int len = coords->length();
  this->_AtomTableIndices->mapHash( [len,&coords] (core::T_sp key, core::T_sp value) {
    int indexX3 = value.unsafe_fixnum()*3;
#ifdef DEBUG_ASSERT
    if ((indexX3%3) || indexX3>=len) SIMPLE_ERROR("indexX3 %d not multiple of 3 or out of range %d", indexX3, len);
#endif
    Atom_sp atom = gc::As_unsafe<Atom_sp>(key);
    Vector3 pos = atom->getPosition();
    (*coords)[indexX3+0] = pos.getX();
    (*coords)[indexX3+1] = pos.getY();
    (*coords)[indexX3+2] = pos.getZ();
  });
}

CL_DOCSTRING("Write coordinates into atoms from an NVector of coordinates");
CL_DEFMETHOD void AtomTable_O::writeAtomCoordinates(NVector_sp coords)
{
  if (coords->length() != this->getNVectorSize()) {
    SIMPLE_ERROR("The passed NVector has length %d but the atom-table needs length %d", coords->length(), this->getNVectorSize() );
  }
  int len = coords->length();
  this->_AtomTableIndices->mapHash( [len,&coords] (core::T_sp key, core::T_sp value) {
    int indexX3 = value.unsafe_fixnum()*3;
#ifdef DEBUG_ASSERT
    if ((indexX3%3) || indexX3>=len) SIMPLE_ERROR("indexX3 %d not multiple of 3 or out of range %d", indexX3, len);
#endif
    Atom_sp atom = gc::As_unsafe<Atom_sp>(key);
    Vector3 pos ((*coords)[indexX3+0],
                 (*coords)[indexX3+1],
                 (*coords)[indexX3+2]);
    atom->setPosition(pos);
  });
}


      
CL_DEFMETHOD size_t AtomTable_O::getCoordinateIndexForAtomAtIndex(size_t index)
{
  if (index<this->_Atoms.size()) {
    EnergyAtom& ea = this->_Atoms[index];
    return this->getCoordinateIndex(ea._SharedAtom);
  }
  SIMPLE_ERROR(("Atom index %d is out of range (0...%d)") , index, this->_Atoms.size());
}

CL_DEFMETHOD core::HashTableEq_sp AtomTable_O::getAtomTableIndices() {
  return this->_AtomTableIndices;
}
CL_DOCSTRING(R"dx(Return the index of ATOM or NIL if it is not found.  The second return value is T if found and NIL if not.)dx");
CL_DEFMETHOD core::T_mv AtomTable_O::getAtomIndexOrNil(Atom_sp atom)
{
  return _AtomTableIndices->gethash(atom);
}

EnergyAtom* AtomTable_O::getEnergyAtomPointer(Atom_sp a)
{
  core::T_mv it = this->_AtomTableIndices->gethash(a);
  core::MultipleValues &values = core::lisp_multipleValues();
  if ( values.second(it.number_of_values()).nilp() ) // it == this->_AtomTableIndices.end() )
  {
    SIMPLE_ERROR(("Could not find atom[%s] in AtomTable") , _rep_(a) );
  }
  return &this->_Atoms[core::clasp_to_fixnum(it)];
}

CL_LISPIFY_NAME("addAtomInfo");
CL_DEFMETHOD     int AtomTable_O::addAtomInfo(Atom_sp atom, units::Quantity_sp charge, units::Quantity_sp mass, int typeIndex, uint atomicNumber )
{
  int coordinateIndex = this->getNumberOfAtoms()*3;
  EnergyAtom ea(atom,coordinateIndex);
  ea._Charge = charge->value_in_unit_asReal(units::_sym_elementaryCharge->symbolValue().as<units::Unit_O>());
  ea._Mass = mass->value_in_unit_asReal(units::_sym_daltons->symbolValue().as<units::Unit_O>());
  ea._TypeIndex = typeIndex;
  ea._AtomicNumber = atomicNumber;
  this->add(ea);
  return coordinateIndex;
}


#ifdef XML_ARCHIVE
void	AtomTable_O::archiveBase(core::ArchiveP node)
{
  if ( node->loading() )
  {
    this->_Atoms.clear();
    core::VectorNodes::iterator	vi;
    for (vi=node->begin_Children(); vi!=node->end_Children(); vi++ )
    {
      EnergyAtom ea(_lisp);
      ea.archive(*vi);
      this->_Atoms.push_back(ea);
      (*vi)->setRecognized(true);
    }
  } else
  {
    core::ArchiveP child;
    gctools::Vec0<EnergyAtom>::iterator	eai;
    for ( eai=this->_Atoms.begin(); eai!=this->_Atoms.end(); eai++ )
    {
      child = node->createChildNode("EnergyAtom");
      eai->archive(child);
    }
  }
}
#endif








CL_DEFMETHOD void	AtomTable_O::dumpTerms(core::HashTable_sp atomTypes)
{
  gctools::Vec0<EnergyAtom>::iterator	eai;
  string				as1,as2,as3,as4;
  core::Symbol_sp				str1, str2, str3, str4;
  int index = 0;
  for ( eai=this->_Atoms.begin(); eai!=this->_Atoms.end(); eai++ ) {
    as1 = atomLabel(eai->atom());
    str1 = eai->atom()->getType(atomTypes);
    core::write_bf_stream(fmt::sprintf("(TERM %d ATOM %-9s %-9s :charge %8.5lf :mass %8.5lf :typeIndex %d)\n"
              , index
              , as1
              , _rep_(str1)
              , eai->_Charge
              , eai->_Mass
                                       , eai->_TypeIndex ));
    index++;
  }
}

CL_DEFMETHOD core::T_sp AtomTable_O::firstSolventMoleculeNSPSOL() const {
  if (this->_firstSolventMoleculeNSPSOL.unboundp()) SIMPLE_ERROR(("The firstSolventMoleculeNSPSOL is unbound"));
  return this->_firstSolventMoleculeNSPSOL;
}
CL_DEFMETHOD bool AtomTable_O::firstSolventMoleculeNSPSOL_boundP() const {
  return !this->_firstSolventMoleculeNSPSOL.unboundp();
}
CL_DEFMETHOD void AtomTable_O::set_firstSolventMoleculeNSPSOL(size_t num) {
  this->_firstSolventMoleculeNSPSOL = core::make_fixnum(num);
}
CL_DEFMETHOD void AtomTable_O::makUnbound_firstSolventMoleculeNSPSOL() {
  this->_firstSolventMoleculeNSPSOL = unbound<core::Fixnum_I>();
}


CL_DEFMETHOD core::T_sp AtomTable_O::finalSoluteResidueIPTRES() const {
  if (this->_finalSoluteResidueIPTRES.unboundp()) SIMPLE_ERROR(("The finalSoluteResidueIPTRES is unbound"));
  return this->_finalSoluteResidueIPTRES;
}
CL_DEFMETHOD bool AtomTable_O::finalSoluteResidueIPTRES_boundP() const {
  return !this->_finalSoluteResidueIPTRES.unboundp();
}
CL_DEFMETHOD void AtomTable_O::set_finalSoluteResidueIPTRES(size_t num) {
  this->_finalSoluteResidueIPTRES = core::make_fixnum(num);
}
CL_DEFMETHOD void AtomTable_O::makUnbound_finalSoluteResidueIPTRES() {
  this->_finalSoluteResidueIPTRES = unbound<core::Fixnum_I>();
}


CL_DEFMETHOD core::T_sp AtomTable_O::totalNumberOfMoleculesNSPM() const {
  if (this->_totalNumberOfMoleculesNSPM.unboundp()) SIMPLE_ERROR(("The totalNumberOfMoleculesNSPM is unbound"));
  return this->_totalNumberOfMoleculesNSPM;
}
CL_DEFMETHOD bool AtomTable_O::totalNumberOfMoleculesNSPM_boundP() const {
  return !this->_totalNumberOfMoleculesNSPM.unboundp();
}
CL_DEFMETHOD void AtomTable_O::set_totalNumberOfMoleculesNSPM(size_t num) {
  this->_totalNumberOfMoleculesNSPM = core::make_fixnum(num);
}
CL_DEFMETHOD void AtomTable_O::makUnbound_totalNumberOfMoleculesNSPM() {
  this->_totalNumberOfMoleculesNSPM = unbound<core::Fixnum_I>();
}






CL_DEFMETHOD core::T_sp AtomTable_O::aggregateName() const {
  if (this->_AggregateName.unboundp()) {
    SIMPLE_ERROR(("The aggregate-name is unbound"));
  }
  return this->_AggregateName;
}

CL_DEFMETHOD bool AtomTable_O::aggregateNameBoundP() const {
  return !this->_AggregateName.unboundp();
}

CL_DEFMETHOD void AtomTable_O::setAggregateName(core::T_sp name) {
  this->_AggregateName = name;
}

CL_DEFMETHOD void AtomTable_O::makUnboundAggregateName() {
  this->_AggregateName = unbound<core::T_O>();
}

SYMBOL_EXPORT_SC_(ChemPkg,bounding_box);
SYMBOL_EXPORT_SC_(ChemPkg,bounding_box_bound_p);
SYMBOL_EXPORT_SC_(ChemPkg,set_bounding_box);
SYMBOL_EXPORT_SC_(ChemPkg,mak_unbound_bounding_box);

CL_DOCSTRING(R"dx(Return the bounding-box for the atom-table.)dx");
CL_LISPIFY_NAME(atom-table-bounding-box);
CL_DEFMETHOD core::T_sp AtomTable_O::boundingBox() const {
  if (this->_BoundingBox.unboundp()) {
    SIMPLE_ERROR(("The bounding-box is unbound"));
  }
  return this->_BoundingBox;
}

CL_DOCSTRING(R"dx(Return T if the bounding-box is bound)dx");
CL_LISPIFY_NAME(atom-table-bounding-box-bound-p);
CL_DEFMETHOD bool AtomTable_O::boundingBoxBoundP() const {
  return !this->_BoundingBox.unboundp();
}

CL_DOCSTRING(R"dx(Set the bounding-box)dx");
CL_LISPIFY_NAME(atom-table-set-bounding-box);
CL_DEFMETHOD void AtomTable_O::setBoundingBox(core::T_sp boundingBox) {
  this->_BoundingBox = boundingBox;
}

CL_DOCSTRING(R"dx(Make the bounding-box unbound)dx");
CL_LISPIFY_NAME(atom-table-mak-unbound-bounding-box);
CL_DEFMETHOD void AtomTable_O::makUnboundBoundingBox() {
  this->_BoundingBox = unbound<core::T_O>();
}

CL_DEFMETHOD core::T_sp AtomTable_O::nonbondForceFieldForAggregate() const {
  if (this->_NonbondForceFieldForAggregate.unboundp()) {
    SIMPLE_ERROR(("The nonbond-forcefield-for-aggregate is unbound"));
  }
  return this->_NonbondForceFieldForAggregate;
}

CL_DEFMETHOD bool AtomTable_O::nonbondForceFieldForAggregateBoundP() const {
  return !this->_NonbondForceFieldForAggregate.unboundp();
}

CL_DEFMETHOD void AtomTable_O::setNonbondForceFieldForAggregate(core::T_sp forceField) {
  this->_NonbondForceFieldForAggregate = forceField;
}

CL_DEFMETHOD void AtomTable_O::makUnboundNonbondForceFieldForAggregate() {
  this->_NonbondForceFieldForAggregate = unbound<FFNonbondDb_O>();
}

CL_DEFMETHOD core::T_sp AtomTable_O::atom_table_residue_pointers() const {
  if (this->_ResiduePointers.unboundp()) {
    SIMPLE_ERROR(("Residue pointers table is not bound"));
  }
  return core::eval::funcall(cl::_sym_copySeq,this->_ResiduePointers);
}

CL_DEFMETHOD core::T_sp AtomTable_O::atom_table_residue_names() const {
  if (this->_ResidueNames.unboundp()) {
    SIMPLE_ERROR(("Residue names table is not bound"));
  }
  return core::eval::funcall(cl::_sym_copySeq,this->_ResidueNames);
}

CL_DEFMETHOD core::T_sp AtomTable_O::atom_table_atoms_per_molecule() const {
  if (this->_AtomsPerMolecule.unboundp()) {
    SIMPLE_ERROR(("atoms per molecule table is not bound"));
  }
  return core::eval::funcall(cl::_sym_copySeq,this->_AtomsPerMolecule);
}

CL_DEFMETHOD core::T_sp AtomTable_O::atom_table_residues() const {
//  printf("%s:%d In :atom_table_residues\n", __FILE__, __LINE__ );
  if (this->_Residues.unboundp()) {
    SIMPLE_ERROR(("residues table is not bound"));
  }
  return core::eval::funcall(cl::_sym_copySeq,this->_Residues);
}

CL_DEFMETHOD core::T_sp AtomTable_O::atom_table_molecules() const {
//  printf("%s:%d In :atom_table_residues\n", __FILE__, __LINE__ );
  if (this->_Molecules.unboundp()) {
    SIMPLE_ERROR(("molecules table is not bound"));
  }
  return core::eval::funcall(cl::_sym_copySeq,this->_Molecules);
}

CL_DEFMETHOD
void AtomTable_O::setAtomFlag(size_t index, size_t flag) {
  if (index>=this->_Atoms.size()) {
    SIMPLE_ERROR(("Index %lu is out of range - must be less than %lu") , index , this->_Atoms.size());
  }
  this->_Atoms[index]._Flag = flag;
}

CL_DEFMETHOD
size_t AtomTable_O::getAtomFlag(size_t index) {
  if (index>=this->_Atoms.size()) {
    SIMPLE_ERROR(("Index %lu is out of range - must be less than %lu") , index , this->_Atoms.size());
  }
  return this->_Atoms[index]._Flag;
}

CL_DEFMETHOD
void AtomTable_O::constructFromMolecule(Molecule_sp mol, core::T_sp nonbondForceField, core::T_sp activeAtoms, core::HashTable_sp atomTypes )
{
  uint idx = this->_Atoms.size();
  uint coordinateIndex = idx*3;
  // Push to the molecules
  this->_Molecules->vectorPushExtend(mol);
  Loop residue_loop;
  residue_loop.loopTopGoal(mol,RESIDUES);
  while (residue_loop.advanceLoopAndProcess()) {
    Residue_sp res = residue_loop.getResidue();
    this->_ResiduePointers->vectorPushExtend(idx);
    this->_ResidueNames->vectorPushExtend(res->getName());
    this->_Residues->vectorPushExtend(res);
    {
      Loop loop;
      Atom_sp a1;
      loop.loopTopGoal(res,ATOMS);
      while ( loop.advanceLoopAndProcess() )
      {
        a1 = loop.getAtom();
        if ( activeAtoms.notnilp() && !inAtomSet(activeAtoms,a1) ) continue;
        if ( a1.isA<VirtualAtom_O>() )
        {
          LOG("Skipping virtual atom[%s]" , _rep_(a1) );
          continue; // skip virtuals
        }
        LOG("Setting atom[%s] in AtomTable[%d]" , _rep_(a1) , idx );
        this->_AtomTableIndices->setf_gethash(a1,core::clasp_make_fixnum(idx));
        EnergyAtom ea(nonbondForceField,a1,coordinateIndex,atomTypes);
        ea._AtomName = a1->getName();
        {
          LOG("Spanning tree for atom: %s\n" , _rep_(a1->getName()));
          SpanningLoop_sp span = SpanningLoop_O::create(a1);
          Atom_sp bonded;
          while ( span->advance() ) {
            bonded = span->getAtom();
            if ( bonded == a1 ) continue;
            int backCount = span->getBackCount(bonded); // ->getBackCount();
            LOG("Looking at atom[%s] at remove[%d]" , _rep_(bonded) , backCount );
			// Once we crawl out 4 bonds we have gone as far as we need
            if ( backCount >= 4 ) {
              LOG("Hit remove of 4 - terminating spanning loop");
              break;
            }
            ASSERT(backCount>0 && backCount<=3);
            LOG("Adding atom at remove %d --> %s\n" , (backCount-1) , _rep_(bonded->getName()));
            ea._AtomsAtRemoveBondAngle14[backCount-1].insert(bonded);
          }
        }
        this->_Atoms.push_back(ea);
        idx++;
        coordinateIndex += 3;
#ifdef DEBUG_ON
        stringstream ss;
        gctools::SmallOrderedSet<Atom_sp>::iterator si;
        for ( int zr = 0; zr<=EnergyAtom::max_remove; zr++ ) {
          ss.str("");
          for ( si = ea._AtomsAtRemoveBondAngle14[zr].begin(); si!=ea._AtomsAtRemoveBondAngle14[zr].end(); si++ ) {
            ss << " " << _rep_((*si));
          }
          LOG("Atoms at remove %d = %s" , zr , ss.str() );
        }
#endif
      }
    }
  }
  this->_ResiduePointers->vectorPushExtend(idx);
  this->_AtomsPerMolecule->vectorPushExtend(idx);
#if 0
  if (idx<3000) {
    printf("%s:%d:%s  mol -> %s idx -> %lu\n", __FILE__, __LINE__, __FUNCTION__, _rep_(mol).c_str(), idx);
  }
#endif
}

/*! Fill excludedAtomIndices with the excluded atom list.
Amber starts counting atoms at 1 so add 1 to every index.
The atomIndex passed is index0.*/
CL_DEFMETHOD size_t AtomTable_O::push_back_excluded_atom_indices_and_sort( core::ComplexVector_int32_t_sp excludedAtomIndices, size_t atomIndex)
{
  size_t start_size = excludedAtomIndices->length();
  EnergyAtom* ea = &(this->_Atoms[atomIndex]);
  uint otherIndex;
  for ( int ri = 0; ri<=EnergyAtom::max_remove; ++ri ) {
    for (auto bi = ea->_AtomsAtRemoveBondAngle14[ri].begin(); bi!=ea->_AtomsAtRemoveBondAngle14[ri].end(); ++bi ) {
      otherIndex = this->_AtomTableIndices->gethash(*bi).unsafe_fixnum();
      // Amber starts counting atom indices from 1 but Clasp starts with 0
      if (otherIndex > atomIndex) {
        excludedAtomIndices->vectorPushExtend(otherIndex);
      }
    }
  }
  size_t end_size = excludedAtomIndices->length();
  if (end_size == start_size ) {
    // Amber rules are that if there are no excluded atoms then put -1 in the
    // excluded atom list - this will be incremented to zero (0) when written out.
    // Oct 2019 -
    //     This is in the Swails description http://ambermd.org/prmtop.pdf
    excludedAtomIndices->vectorPushExtend(-1); // Use -1 here because it will be incremented by one when writing out and will be written as 0
    ++end_size;
  }
  // sort the indices in increasing order
  sort::quickSortMemory((int32_t*)excludedAtomIndices->rowMajorAddressOfElement_(0),start_size,end_size);
  if (chem__verbose(1)) {
    core::write_bf_stream(fmt::sprintf("%d: " , atomIndex));
    for ( size_t ii = start_size; ii<end_size; ++ii ) {
      core::write_bf_stream(fmt::sprintf("%d " , (*excludedAtomIndices)[ii]));
    }
    core::write_bf_stream(fmt::sprintf("\n"));
  }
  return (end_size - start_size);
}

SYMBOL_EXPORT_SC_(ClPkg,copy_seq);

/*! Calculate the AMBER excluded atom list and return two vectors, one containing the number of
excluded atoms for each atom and the second containing the sorted excluded atom list */
CL_DEFMETHOD core::T_mv AtomTable_O::calculate_excluded_atom_list()
{
//  printf("%s:%d In calculate_excludec_atom_list\n", __FILE__, __LINE__ );

  core::ComplexVector_int32_t_sp number_excluded_atoms = core::ComplexVector_int32_t_O::make_vector(32,0,core::make_fixnum(0),nil<core::T_O>(),false,core::make_fixnum(0));
  core::ComplexVector_int32_t_sp excluded_atoms_list = core::ComplexVector_int32_t_O::make_vector(32,0,core::make_fixnum(0),nil<core::T_O>(),false,core::make_fixnum(0));
  size_t num_atoms = this->getNumberOfAtoms();
  for ( size_t i1=0; i1<num_atoms; ++i1) {
    size_t num = this->push_back_excluded_atom_indices_and_sort(excluded_atoms_list,i1);
    number_excluded_atoms->vectorPushExtend(num);
  }
  core::T_sp t_number_excluded_atoms = core::eval::funcall(cl::_sym_copy_seq,number_excluded_atoms);
  core::T_sp t_excluded_atoms_list = core::eval::funcall(cl::_sym_copy_seq,excluded_atoms_list);
  return Values(t_number_excluded_atoms, t_excluded_atoms_list);
}

CL_DEFMETHOD core::Symbol_sp AtomTable_O::elt_atom_type(int index,core::HashTable_sp atomTypes) {
  return this->_Atoms[index]._SharedAtom->getType(atomTypes);
};
SYMBOL_EXPORT_SC_(KeywordPkg,atom_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,atom_name_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,atom_type_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,charge_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,mass_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,atomic_number_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,residues);
SYMBOL_EXPORT_SC_(KeywordPkg,residue_pointers);
SYMBOL_EXPORT_SC_(KeywordPkg,residue_names);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms_per_molecule);

CL_DEFMETHOD void  AtomTable_O::fill_atom_table_from_vectors(core::List_sp vectors)
{
  core::Array_sp atoms_vec = (safe_alist_lookup<core::Array_sp>(vectors,kw::_sym_atom_vector));
  core::Array_sp atom_name_vec = (safe_alist_lookup<core::Array_sp>(vectors,kw::_sym_atom_name_vector));
  core::Array_sp atom_type_vec = (safe_alist_lookup<core::Array_sp>(vectors,kw::_sym_atom_type_vector));
  core::Array_sp charge_vec = (safe_alist_lookup<core::Array_sp>(vectors,kw::_sym_charge_vector));
  core::Array_sp mass_vec = (safe_alist_lookup<core::Array_sp>(vectors,kw::_sym_mass_vector));
  core::Array_sp atomic_number_vec = (safe_alist_lookup<core::Array_sp>(vectors,kw::_sym_atomic_number_vector));
//  this->_Residues = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_residues));
//  this->_ResidueNames = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_residue_names));
  this->_Atoms.resize(atom_name_vec->length());
  this->_AtomTableIndices->clrhash();
  for (size_t i = 0, iEnd(atom_name_vec->length()); i<iEnd ;++i)
  {
//    printf("%s:%d About to set _AtomName with %s\n", __FILE__, __LINE__, _rep_(atom_name_vec->rowMajorAref(i)).c_str());
    Atom_sp atom = gc::As<Atom_sp>(atoms_vec->rowMajorAref(i));
    this->_AtomTableIndices->setf_gethash(atom,core::clasp_make_fixnum(i));
    core::T_sp type = atom_type_vec->rowMajorAref(i);
//    printf("%s:%d  type -> %s\n", __FILE__, __LINE__, _rep_(type).c_str());
//    atom->setAtomType(atom_type_vec->rowMajorAref(i));
    double charge  =  translate::from_object<double>(charge_vec->rowMajorAref(i))._v;   // charge-vector
    this->_Atoms[i]._Charge       =  charge;
    atom->setCharge(charge);
    this->_Atoms[i]._SharedAtom     =  atom;
    this->_Atoms[i]._AtomName     =  gc::As<core::Symbol_sp>(atom_name_vec->rowMajorAref(i));  // atom-name-vector
    //    The _TypeIndex is going to go away once we have ensured that the new Common Lisp code
    // that calculates nonbond terms works.
//    this->_Atoms[i]._TypeIndex    =  forceField->getNonbondDb()->findTypeIndex((*atom_type_vec)[i]);
    this->_Atoms[i]._TypeIndex = UNDEF_UINT;
    this->_Atoms[i]._Mass         =  translate::from_object<double>(mass_vec->rowMajorAref(i))._v;                // masses
    this->_Atoms[i]._AtomicNumber =  translate::from_object<int>(atomic_number_vec->rowMajorAref(i))._v;       // vec
  }
}


DOCGROUP(cando);
CL_DEFUN core::List_sp chem__atoms_at_remove(AtomTable_sp table, size_t index, size_t remove ) {
  if (index<table->_Atoms.size()) {
    EnergyAtom* ea = &(table->_Atoms[index]);
    if (remove <4) {
      core::T_sp result = nil<core::T_O>();
      for (auto bi = ea->_AtomsAtRemoveBondAngle14[remove].begin(); bi!=ea->_AtomsAtRemoveBondAngle14[remove].end(); ++bi ) {
        result = core::Cons_O::create(*bi,result);
      }
      return result;
    }
    SIMPLE_ERROR(("The remove is out of bounds - %d must be <=3") , remove);
  }
  SIMPLE_ERROR(("Atom index %d is out of bounds - must be < %d") , index , table->_Atoms.size());
}



};
