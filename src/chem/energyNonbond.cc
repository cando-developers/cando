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
#define	DEBUG_LEVEL_FULL

#include <clasp/core/foundation.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/array.h>
#include <clasp/core/evaluator.h>
#include <cando/geom/color.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/nativeVector.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/matter.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/loop.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
//#define LOG_ENERGY core::write_bf_stream

namespace chem
{

core::List_sp EnergyNonbond::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,da),core::clasp_make_double_float(this->term.dA)),
                                  core::Cons_O::create(INTERN_(kw,dc),core::clasp_make_double_float(this->term.dC)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.I2)),
                                  core::Cons_O::create(INTERN_(kw,charge1), core::clasp_make_double_float(this->_Charge1)),
                                  core::Cons_O::create(INTERN_(kw,charge2), core::clasp_make_double_float(this->_Charge2))
                                  //core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                  //core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2)
                                  );
}

 void EnergyNonbond::decode(core::List_sp alist) {
  SIMPLE_ERROR(BF("Implement decode of EnergyNonbond"));
}

SYMBOL_EXPORT_SC_(ChemPkg,find_type);
//
// Return true if we could fill the energyNonbond term
// otherwise false usually if we don't recognize one of the atom types like DU
//
bool	EnergyNonbond::defineFrom(core::T_sp	forceField,
                                  bool		is14,
                                  EnergyAtom	*iea1,
                                  EnergyAtom	*iea2,
                                  EnergyNonbond_sp energyNonbond)
{_OF();
  double				epsilonij;
  double				vdwScale;
  double				electrostaticScale;
  this->_Is14 = is14;
  this->_Atom1 = iea1->atom();
  this->_Atom2 = iea2->atom();
  core::Symbol_sp t1 = iea1->atom()->getType();
  core::Symbol_sp t2 = iea2->atom()->getType();
  LOG(BF("Defining nonbond between types: %s - %s") % _rep_(t1) % _rep_(t2));
  ASSERT(forceField&&forceField.notnilp());
  LOG(BF("forceField @%p   .notnilp()->%d") % forceField.raw_() % forceField.notnilp());
  core::T_sp tffNonbond1 = core::eval::funcall(_sym_find_type,forceField,t1);
  core::T_sp tffNonbond2 = core::eval::funcall(_sym_find_type,forceField,t2);
  ANN(tffNonbond1);
  if ( tffNonbond1.nilp() )
  {
//     	SIMPLE_ERROR(BF("Unknown force field type(",iea1->_Atom->getType().c_str(),") for non-bonded interaction"));
    return false;
  }
  ANN(tffNonbond2);
  if ( tffNonbond2.nilp() )
  {
//     	SIMPLE_ERROR(BF("Unknown force field type(",iea2->_Atom->getType().c_str(),") for non-bonded interaction"));
    return false;
  }
  if ( is14 ) {
    vdwScale = 1.0/2.0*energyNonbond->getVdwScale();
    electrostaticScale = energyNonbond->getElectrostaticScale()*1.0/1.2;
  } else
  {
    vdwScale = energyNonbond->getVdwScale();
    electrostaticScale = energyNonbond->getElectrostaticScale();
  }
  LOG(BF( "vdwScale = %lf")% (double)(vdwScale) );
  LOG(BF( "electrostaticScale = %lf")% (double)(electrostaticScale) );
  LOG(BF( " is14=%d")% is14 );
  FFNonbond_sp ffNonbond1 = gc::As<FFNonbond_sp>(tffNonbond1);
  FFNonbond_sp ffNonbond2 = gc::As<FFNonbond_sp>(tffNonbond2);
  {_BLOCK_TRACE("Calculating nonbond parameters");
    this->_RStar = ffNonbond1->getRadius_Angstroms()+ffNonbond2->getRadius_Angstroms();
    epsilonij = sqrt(ffNonbond1->getEpsilon_kCal()*ffNonbond2->getEpsilon_kCal());
    this->_A = epsilonij*pow(this->_RStar,12.0);
    this->_C = 2.0*epsilonij*pow(this->_RStar,6.0);
    this->term.dA = this->_A*vdwScale;
    this->term.dC = this->_C*vdwScale;
  }
  {_BLOCK_TRACE("Calculating electrostatic parameters");
    electrostaticScale *= ELECTROSTATIC_MODIFIER;
    this->_Charge1 = iea1->atom()->getCharge();
    this->_Charge2 = iea2->atom()->getCharge();
    this->term.dQ1Q2 = electrostaticScale*(this->_Charge1*this->_Charge2)/energyNonbond->getDielectricConstant();
    LOG(BF( "Calc dQ1Q2 electrostaticScale= %lf")% (double)(electrostaticScale));
    LOG(BF( "Calc dQ1Q2 Dielectric constant = %lf")% (double)(energyNonbond->getDielectricConstant()));
    LOG(BF( "Calc dQ1Q2 Charge1 = %lf")% (double)(this->_Charge1));
    LOG(BF( "Calc dQ1Q2 Charge2 = %lf")% (double)(this->_Charge2));
    LOG(BF( "dQ1Q2 = %lf")% (double)(this->term.dQ1Q2));
  }
  this->term.I1 = iea1->coordinateIndexTimes3();
  this->term.I2 = iea2->coordinateIndexTimes3();
  return true;
}

double	EnergyNonbond::getDistance()
{
  Vector3	pos1, pos2;
  pos1 = this->_Atom1->getPosition();
  pos2 = this->_Atom2->getPosition();
  return geom::calculateDistance(pos1,pos2);
}


#if 0
adapt::QDomNode_sp	EnergyNonbond::asXml()
{
  adapt::QDomNode_sp	node;
  Vector3	vdiff;

  node = adapt::QDomNode_O::create(lisp,"EnergyNonbond");
  node->addAttributeString("atom1Name",this->_Atom1->getName());
  node->addAttributeString("atom2Name",this->_Atom2->getName());
  node->addAttributeInt("I1",this->term.I1);
  node->addAttributeInt("I2",this->term.I2);
  node->addAttributeBool("is14",this->_Is14);
  node->addAttributeString("atom1Type",this->_Atom1->getType());
  node->addAttributeString("atom2Type",this->_Atom2->getType());
  node->addAttributeDoubleScientific("RStar",this->_RStar);
  node->addAttributeDoubleScientific("A",this->_A);
  node->addAttributeDoubleScientific("C",this->_C);
  node->addAttributeDoubleScientific("Charge1",this->_Charge1);
  node->addAttributeDoubleScientific("Charge2",this->_Charge2);
//    vdiff = this->_Atom1->_Atom->getPosition() - this->_Atom2->_Atom->getPosition();
//    diff = vdiff.length();
//    node->addAttributeDouble("_r",diff,5,2);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
  adapt::QDomNode_sp xml = adapt::QDomNode_O::create(lisp,"Evaluated");
  xml->addAttributeBool("calcForce",this->_calcForce );
  xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
  xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include <_Nonbond_debugEvalXml.cc>
  node->addChild(xml);
#endif
  node->addAttributeDoubleScientific("dA",this->term.dA);
  node->addAttributeDoubleScientific("dC",this->term.dC);
  node->addAttributeDoubleScientific("dQ1Q2",this->term.dQ1Q2);
  return node;
}

void	EnergyNonbond::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
                                                  AtomTable_sp at)
{
  this->term.dA = xml->getAttributeDouble("dA");
  this->term.dC = xml->getAttributeDouble("dC");
  this->term.dQ1Q2 = xml->getAttributeDouble("dQ1Q2");
  this->_RStar = xml->getAttributeDouble("RStar");
  this->_A = xml->getAttributeDouble("A");
  this->_C = xml->getAttributeDouble("C");
  this->_Charge1 = xml->getAttributeDouble("Charge1");
  this->_Charge2 = xml->getAttributeDouble("Charge2");
  this->term.I1 = xml->getAttributeInt("I1");
  this->term.I2 = xml->getAttributeInt("I2");
  this->_Is14 = xml->getAttributeBool("is14");
  this->_Atom1 = at->findEnergyAtomWithCoordinateIndex(this->term.I1)->atom();
  this->_Atom2 = at->findEnergyAtomWithCoordinateIndex(this->term.I2)->atom();
}
#endif



double	_evaluateEnergyOnly_Nonbond(double x1, double y1, double z1,
                                    double x2, double y2, double z2,
                                    double dA, double dC, double dQ1Q2)
{
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {}
#undef	NONBOND_FORCE_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE(i,o,v) {}
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(x) {}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(x) {}

#undef	NONBOND_CALC_FORCE	// Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>

  return Energy;
}


void	EnergyNonbond_O::zeroEnergy()
{
  this->Base::zeroEnergy();
  this->_EnergyElectrostatic = 0.0;
  this->_EnergyVdw = 0.0;
}

double	EnergyNonbond_O::getEnergy()
{
  double	e;
  e = this->getVdwEnergy();
  e += this->getElectrostaticEnergy();
  return e;
}


void	EnergyNonbond_O::dumpTerms()
{
  gctools::Vec0<EnergyNonbond>::iterator	eni;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4;
  for ( eni=this->_Terms.begin(); eni!=this->_Terms.end(); eni++ )
  {
    as1 = atomLabel(eni->_Atom1);
    as2 = atomLabel(eni->_Atom2);
    if ( as1 < as2 )
    {
      str1 = as1;
      str2 = as2;
    } else {
      str2 = as1;
      str1 = as2;
    }
    if ( eni->_Is14 ) {
      _lisp->print(BF("TERM 4CALC14 %-9s - %-9s") %	str1 % str2 );
    }else{
      _lisp->print(BF("TERM 5NONBOND %-9s - %-9s") %	str1 % str2 );
    }
  }

}


#ifdef RENDER
int EnergyNonbond_O::countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn, core::Lisp_sp lisp)
{
  gctools::Vec0<EnergyNonbond>::iterator	eni;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4;
  int				overlapCount = 0;
  int				ia1, ia2;
  double				x1, y1, z1;
  double				x2, y2, z2;
  double				dx,dy,dz;
  double				cutoff, distSquared;
  bool				render;
  Vector3				v1,v2;
  geom::Color_sp				color;
  geom::GrLines_sp			lines;
  lines = geom::GrLines_O::create();
  render = false;
  ANN(displayIn);
  if ( displayIn.notnilp() )
  {
    render = true;
    displayIn->clear();
    color = geom::Color_O::systemColor(kw::_sym_yellow);
  }
  for ( eni=this->_Terms.begin();
        eni!=this->_Terms.end(); eni++ )
  {
    ia1 = eni->term.I1;
    ia2 = eni->term.I2;
    x1 = pos->element(ia1+0);
    y1 = pos->element(ia1+1);
    z1 = pos->element(ia1+2);
    x2 = pos->element(ia2+0);
    y2 = pos->element(ia2+1);
    z2 = pos->element(ia2+2);
    cutoff = eni->_RStar*scaleSumOfVdwRadii;
    cutoff = cutoff*cutoff;
    dx = x1-x2;
    dy = y1-y2;
    dz = z1-z2;
    distSquared = dx*dx+dy*dy+dz*dz;
    if ( distSquared<cutoff )
    {
      LOG(BF("Found a close contact with distance = %lf and cutoff %lf") % sqrt(distSquared) % sqrt(cutoff)  );
      LOG(BF("Atom1 = %s") % eni->_Atom1->description()  );
      LOG(BF("Atom2 = %s") % eni->_Atom2->description()  );
      if ( render )
      {
        v1.set(x1,y1,z1);
        v2.set(x2,y2,z2);
        lines->appendLineWithVertices(v1,color,v2,color);
      }
      overlapCount++;
    }
  }
  displayIn->append(lines);
  return overlapCount;
}
#endif


void	EnergyNonbond_O::setupHessianPreconditioner(
                                                    NVector_sp nvPosition,
                                                    AbstractLargeSquareMatrix_sp m )
{
  SIMPLE_ERROR(BF("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!!"));
}

double	EnergyNonbond_O::evaluateAll(NVector_sp 	pos,
                                     bool 		calcForce,
                                     gc::Nilable<NVector_sp> 	force,
                                     bool		calcDiagonalHessian,
                                     bool		calcOffDiagonalHessian,
                                     gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                     gc::Nilable<NVector_sp>	hdvec, 
                                     gc::Nilable<NVector_sp> 	dvec )
{
//  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  if (this->_UsesExcludedAtoms) {
    // Evaluate the nonbonds using the excluded atom list
    this->evaluateUsingExcludedAtoms(pos,calcForce,force,calcDiagonalHessian,
                                     calcOffDiagonalHessian,hessian,hdvec,dvec);
    // Evaluate the 1-4 terms
    this->evaluateTerms(pos,calcForce,force,calcDiagonalHessian,
                             calcOffDiagonalHessian,hessian,hdvec,dvec);
  } else {
    // Evaluate everything using terms
    this->evaluateTerms(pos,calcForce,force,calcDiagonalHessian,
                             calcOffDiagonalHessian,hessian,hdvec,dvec);
  }
  this->_TotalEnergy = this->_EnergyVdw+this->_EnergyElectrostatic;
  return this->_TotalEnergy;
}
    
    

__attribute__((optnone)) void	EnergyNonbond_O::evaluateTerms(NVector_sp 	pos,
                                     bool 		calcForce,
                                     gc::Nilable<NVector_sp> 	force,
                                     bool		calcDiagonalHessian,
                                     bool		calcOffDiagonalHessian,
                                     gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                     gc::Nilable<NVector_sp>	hdvec, 
                                     gc::Nilable<NVector_sp> 	dvec )
{
//  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{x=nbi->term.x;}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);}

#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
  if ( this->isEnabled() ) 
  {
    gctools::Vec0<EnergyNonbond>::iterator firstElement = this->_Terms.begin();
    int nonBondTerms = this->_Terms.size();
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
    {
      LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>

#pragma clang diagnostic pop

      double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
      int	I1, I2,i;
      gctools::Vec0<EnergyNonbond>::iterator nbi;
      for ( i=0;i<nonBondTerms; i++ )
      {
        nbi = firstElement+i;
//        printf("%s:%d:%s i %d \n", __FILE__, __LINE__, __FUNCTION__,i );

#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
        if ( this->_Debug_NumberOfNonbondTermsToCalculate > 0 ) {
  printf("%s:%d:%s \n", __FILE__, __LINE__, __FUNCTION__ );
          if ( i>= this->_Debug_NumberOfNonbondTermsToCalculate ) {
            break;
          }
        }
#endif
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>

#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
        nbi->_calcForce = calcForce;
        nbi->_calcDiagonalHessian = calcDiagonalHessian;
        nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET

#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/chem/energy_functions/_Nonbond_debugEvalSet.cc>
#endif //]

#ifdef DEBUG_NONBOND_TERM
        if ( this->_DebugEnergy ) {
          std::string key;
          std::string key1 = nbi->_Atom1->getName()->symbolNameAsString();
          std::string key2 = nbi->_Atom2->getName()->symbolNameAsString();
          if ( key1<key2 ) {
            key = key1+"-"+key2;
          } else {
            key = key2+"-"+key1;
          }
          if (nbi->_Is14 ) {
            LOG_ENERGY(BF( "MEISTER nonbond %s 1-4 args cando\n")% key );
          } else {
            LOG_ENERGY(BF( "MEISTER nonbond %s args cando\n")% key );
          }
          LOG_ENERGY(BF( "MEISTER nonbond %s dA %5.3lf\n")% key % (nbi->term.dA) );
          LOG_ENERGY(BF( "MEISTER nonbond %s dC %5.3lf\n")% key % (nbi->term.dC) );
          LOG_ENERGY(BF( "MEISTER nonbond %s dQ1Q2 %5.3lf\n")% key % (nbi->term.dQ1Q2) );
          LOG_ENERGY(BF( "MEISTER nonbond %s x1 %5.3lf %d\n")% key % x1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s y1 %5.3lf %d\n")% key % y1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s z1 %5.3lf %d\n")% key % z1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s x2 %5.3lf %d\n")% key % x2 % (I2/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s y2 %5.3lf %d\n")% key % y2 % (I2/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s z2 %5.3lf %d\n")% key % z2 % (I2/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s results\n")% key );
          LOG_ENERGY(BF( "MEISTER nonbond %s evdw %lf\n")% key % Evdw);
          LOG_ENERGY(BF( "MEISTER nonbond %s eeel %lf\n")% key % Eeel);
          LOG_ENERGY(BF( "MEISTER nonbond %s Enonbond(evdw+eeel) %lf\n")% key % (Evdw+Eeel) );
          if ( calcForce ) {
            LOG_ENERGY(BF( "MEISTER nonbond %s fx1 %lf %d\n")% key % fx1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s fy1 %lf %d\n")% key % fy1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s fz1 %lf %d\n")% key % fz1 % (I1/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s fx2 %lf %d\n")% key % fx2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s fy2 %lf %d\n")% key % fy2 % (I2/3+1) );
            LOG_ENERGY(BF( "MEISTER nonbond %s fz2 %lf %d\n")% key % fz2 % (I2/3+1) );
          }
          LOG_ENERGY(BF( "MEISTER nonbond %s stop\n")% key );
        }
#endif
      }
    }
  } else {
    LOG(BF("Nonbond component is not enabled"));
  }
//  printf( "Nonbond energy vdw(%lf) electrostatic(%lf)\n", (double)this->_EnergyVdw, this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF("Nonbond energy }\n"));
}



__attribute__((optnone)) void	EnergyNonbond_O::evaluateUsingExcludedAtoms(NVector_sp 	pos,
                                                    bool 		calcForce,
                                                    gc::Nilable<NVector_sp> 	force,
                                                    bool		calcDiagonalHessian,
                                                    bool		calcOffDiagonalHessian,
                                                    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                    gc::Nilable<NVector_sp>	hdvec, 
                                                    gc::Nilable<NVector_sp> 	dvec )
{
//  printf("%s:%d In evaluateUsingExcludedAtoms starting this->_DebugEnergy -> %d\n", __FILE__, __LINE__, this->_DebugEnergy );
  if (!this->_iac_vec) {
    SIMPLE_ERROR(BF("The nonbonded excluded atoms parameters have not been set up"));
  }
  core::SimpleVector_int32_t_sp numberOfExcludedAtoms = this->_NumberOfExcludedAtomIndices;
  core::SimpleVector_int32_t_sp excludedAtomIndices = this->_ExcludedAtomIndices;
  double vdwScale = this->getVdwScale();
  double electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
//  printf("%s:%d electrostaticcharge %lf\n", __FILE__, __LINE__, electrostaticScale );
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
//  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  if ( !this->isEnabled() ) return;
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
  LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
  // printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2,dA_old,dC_old,dQ1Q2_old;
  int	I1, I2;
  int i = 0;
  int maxIndex = pos->length()/3;
//  int maxIndex = this->_AtomTable->getNumberOfAtoms();
  int excludedAtomIndex = 0;
  int nlocaltype = 0;
  // Find the max local type
  for (i=0; i<this->_iac_vec->length(); ++i){
    if (nlocaltype < (*this->_iac_vec)[i]){
      nlocaltype = (*this->_iac_vec)[i];
        }
  }

  for ( int index1 = 0, index1_end(maxIndex-1); index1 <index1_end; ++index1 ) {
    LOG(BF("%s ====== top of outer loop - index1 = %d\n") % __FUNCTION__ % index1 );
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    LOG(BF("Read numberOfExcludedAtomsRemaining from numberOfExcludedAtoms[%d]= %d\n") % index1 % numberOfExcludedAtomsRemaining);
          // Skip 0 in excluded atom list that amber requires
    if (numberOfExcludedAtomsRemaining<0) {
      ++excludedAtomIndex;
    }
    double charge11 = (*this->_charge_vector)[index1];
    double electrostatic_scaled_charge11 = charge11*electrostaticScale;
    for ( int index2 = index1+1, index2_end(maxIndex); index2 < index2_end; ++index2 ) {
      LOG(BF("    --- top of inner loop   numberOfExcludedAtomsRemaining -> %d    index2 -> %d\n") % numberOfExcludedAtomsRemaining % index2 );
      if (numberOfExcludedAtomsRemaining>0) {
        LOG(BF("    excludedAtomIndices[%d] -> %d  index2 -> %d\n") % excludedAtomIndex % (*excludedAtomIndices)[excludedAtomIndex] % index2 );
        if ((*excludedAtomIndices)[excludedAtomIndex] == index2) {
          LOG(BF("    Excluding atom %d\n") % index2);
          ++excludedAtomIndex;
          --numberOfExcludedAtomsRemaining;
          continue;
        }
      }
      //   cn1i    (*this->_cn1_vec)[cn1i]
      int localindex1 = (*this->_iac_vec)[index1];
      int localindex2 = (*this->_iac_vec)[index2];
      dA = (*this->_cn1_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
      dC = (*this->_cn2_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
//      printf("%s:%d localindex1 %d and localindex2 %d\n", __FILE__, __LINE__, localindex1, localindex2);
//      printf("%s:%d dA     %lf and dC     %lf\n", __FILE__, __LINE__, dA, dC);
      double charge22 = (*this->_charge_vector)[index2];
//      dQ1Q2 = electrostatic_scaled_charge11*charge22;
      dQ1Q2 = charge11*charge22;
//      printf("%s:%d charge1     %lf and charge2     %lf\n", __FILE__, __LINE__, charge11, charge22);
//      printf("%s:%d electrostaticScale     %lf and dQ1Q2     %lf\n", __FILE__, __LINE__, electrostaticScale, dQ1Q2);
#ifdef DEBUG_NONBOND_TERM
      if ( this->_DebugEnergy ) {
        LOG_ENERGY(BF( "nonbond localindex1 %d\n")% localindex1 );
        LOG_ENERGY(BF( "nonbond localindex2 %d\n")% localindex2 );
        LOG_ENERGY(BF( "nonbond dA %5.3lf\n")% dA );
        LOG_ENERGY(BF( "nonbond dC %5.3lf\n")% dC );
        LOG_ENERGY(BF( "nonbond dQ1Q2 %5.3lf\n")% dQ1Q2 );
      }
#endif
      ////////////////////////////////////////////////////////////
      //
      // To here
      //
      ////////////////////////////////////////////////////////////
      I1 = index1*3; 
      I2 = index2*3; 
#ifdef DEBUG_NONBOND_TERM
      if ( this->_DebugEnergy ) {
        LOG_ENERGY(BF( "nonbond I1 %d\n")% I1 );
        LOG_ENERGY(BF( "nonbond I2 %d\n")% I2 );
      }
#endif
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfNonbondTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfNonbondTermsToCalculate ) {
          break;
        }
      }
#endif
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
      nbi->_calcForce = calcForce;
      nbi->_calcDiagonalHessian = calcDiagonalHessian;
      nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/chem/energy_functions/_Nonbond_debugEvalSet.cc>
#endif //]
#ifdef DEBUG_NONBOND_TERM
      if ( this->_DebugEnergy ) {
        std::string key;
        std::string atom1Name = gc::As<core::Symbol_sp>((*this->_atom_name_vector)[index1])->symbolNameAsString();
        std::string atom2Name = gc::As<core::Symbol_sp>((*this->_atom_name_vector)[index2])->symbolNameAsString();
        if ( atom1Name < atom2Name ) {
          key = atom1Name+"-"+atom2Name;
        } else {
          key = atom2Name+"-"+atom1Name;
        }
        LOG_ENERGY(BF( "MEISTER nonbond %s args cando\n")% key );
        LOG_ENERGY(BF( "MEISTER nonbond %s dA %5.3lf\n")% key % dA );
        LOG_ENERGY(BF( "MEISTER nonbond %s dC %5.3lf\n")% key % dC );
        LOG_ENERGY(BF( "MEISTER nonbond %s dQ1Q2 %5.3lf\n")% key % dQ1Q2 );
        LOG_ENERGY(BF( "MEISTER nonbond %s x1 %5.3lf %d\n")% key % x1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s y1 %5.3lf %d\n")% key % y1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s z1 %5.3lf %d\n")% key % z1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s x2 %5.3lf %d\n")% key % x2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s y2 %5.3lf %d\n")% key % y2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s z2 %5.3lf %d\n")% key % z2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s results\n")% key );
        LOG_ENERGY(BF( "MEISTER nonbond %s evdw %lf\n")% key % Evdw);
        LOG_ENERGY(BF( "MEISTER nonbond %s eeel %lf\n")% key % Eeel);
        LOG_ENERGY(BF( "MEISTER nonbond %s Enonbond(evdw+eeel) %lf\n")% key % (Evdw+Eeel) );
        if ( calcForce ) {
          LOG_ENERGY(BF( "MEISTER nonbond %s fx1 %lf %d\n")% key % fx1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s fy1 %lf %d\n")% key % fy1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s fz1 %lf %d\n")% key % fz1 % (I1/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s fx2 %lf %d\n")% key % fx2 % (I2/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s fy2 %lf %d\n")% key % fy2 % (I2/3+1) );
          LOG_ENERGY(BF( "MEISTER nonbond %s fz2 %lf %d\n")% key % fz2 % (I2/3+1) );
        }
        LOG_ENERGY(BF( "MEISTER nonbond %s stop\n")% key );
      }
#endif
    }
  }
//  printf( "Nonbond energy vdw(%lf) electrostatic(%lf)\n", (double)this->_EnergyVdw,  this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy }\n"));
}




CL_DEFMETHOD void EnergyNonbond_O::expandExcludedAtomsToTerms()
{
//  printf("%s:%d In evaluateUsingExcludedAtoms starting this->_DebugEnergy -> %d\n", __FILE__, __LINE__, this->_DebugEnergy );
  if (!this->_iac_vec) {
    SIMPLE_ERROR(BF("The nonbonded excluded atoms parameters have not been set up"));
  }
  core::SimpleVector_int32_t_sp numberOfExcludedAtoms = this->_NumberOfExcludedAtomIndices;
  core::SimpleVector_int32_t_sp excludedAtomIndices = this->_ExcludedAtomIndices;
  double vdwScale = this->getVdwScale();
  double electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
  LOG(BF("Nonbond component is enabled") );
  // printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2,dA_old,dC_old,dQ1Q2_old;
  int	I1, I2;
  int i = 0;
//  int maxIndex = pos->length()/3;
  int maxIndex = this->_AtomTable->getNumberOfAtoms();
  printf( "maxindex %d\n", maxIndex);
  int nlocaltype = 0;
  // Find the max local type
  for (i=0; i<this->_iac_vec->length(); ++i){
    if (nlocaltype < (*this->_iac_vec)[i]){
      nlocaltype = (*this->_iac_vec)[i];
    }
  }
  printf("%s:%d:%s    values -> %lu\n", __FILE__, __LINE__, __FUNCTION__, this->_Terms.size());

  // Count the number of nonbond cross terms that we are going to generate
  size_t excludedAtomIndex = 0;
  size_t count = 0;
  for ( int index1 = 0, index1_end(maxIndex-1); index1 <index1_end; ++index1 ) {
    LOG(BF("%s ====== top of outer loop - index1 = %d\n") % __FUNCTION__ % index1 );
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    LOG(BF("Read numberOfExcludedAtomsRemaining from numberOfExcludedAtoms[%d]= %d\n") % index1 % numberOfExcludedAtomsRemaining);
          // Skip 0 in excluded atom list that amber requires
    if (numberOfExcludedAtomsRemaining<0) {
      ++excludedAtomIndex;
    }
    for ( int index2 = index1+1, index2_end(maxIndex); index2 < index2_end; ++index2 ) {
      LOG(BF("    --- top of inner loop   numberOfExcludedAtomsRemaining -> %d    index2 -> %d\n") % numberOfExcludedAtomsRemaining % index2 );
      if (numberOfExcludedAtomsRemaining>0) {
        LOG(BF("    excludedAtomIndices[%d] -> %d  index2 -> %d\n") % excludedAtomIndex % (*excludedAtomIndices)[excludedAtomIndex] % index2 );
        if ((*excludedAtomIndices)[excludedAtomIndex] == index2) {
          LOG(BF("    Excluding atom %d\n") % index2);
          ++excludedAtomIndex;
          --numberOfExcludedAtomsRemaining;
          continue;
        }
      }
      ++count;
    }
  }
  this->_Terms.resize(count);
  printf("%s:%d:%s    values -> %lu\n", __FILE__, __LINE__, __FUNCTION__, this->_Terms.size());

  // Now fill in the terms
  excludedAtomIndex = 0;
  size_t termIndex = 0;
  for ( int index1 = 0, index1_end(maxIndex-1); index1 <index1_end; ++index1 ) {
    LOG(BF("%s ====== top of outer loop - index1 = %d\n") % __FUNCTION__ % index1 );
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    LOG(BF("Read numberOfExcludedAtomsRemaining from numberOfExcludedAtoms[%d]= %d\n") % index1 % numberOfExcludedAtomsRemaining);
          // Skip 0 in excluded atom list that amber requires
    if (numberOfExcludedAtomsRemaining<0) {
      ++excludedAtomIndex;
    }
    double charge11 = (*this->_charge_vector)[index1];
    double electrostatic_scaled_charge11 = charge11*electrostaticScale;
    for ( int index2 = index1+1, index2_end(maxIndex); index2 < index2_end; ++index2 ) {
      LOG(BF("    --- top of inner loop   numberOfExcludedAtomsRemaining -> %d    index2 -> %d\n") % numberOfExcludedAtomsRemaining % index2 );
      if (numberOfExcludedAtomsRemaining>0) {
        LOG(BF("    excludedAtomIndices[%d] -> %d  index2 -> %d\n") % excludedAtomIndex % (*excludedAtomIndices)[excludedAtomIndex] % index2 );
        if ((*excludedAtomIndices)[excludedAtomIndex] == index2) {
          LOG(BF("    Excluding atom %d\n") % index2);
          ++excludedAtomIndex;
          --numberOfExcludedAtomsRemaining;
          continue;
        }
      }
      //   cn1i    (*this->_cn1_vec)[cn1i]
      int localindex1 = (*this->_iac_vec)[index1];
      int localindex2 = (*this->_iac_vec)[index2];
      dA = (*this->_cn1_vec)[(*this->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
      dC = (*this->_cn2_vec)[(*this ->_ico_vec)[nlocaltype*(localindex1-1)+localindex2-1]-1];
//      printf("%s:%d localindex1 %d and localindex2 %d\n", __FILE__, __LINE__, localindex1, localindex2);
//      printf("%s:%d dA     %lf and dC     %lf\n", __FILE__, __LINE__, dA, dC);
      double charge22 = (*this->_charge_vector)[index2];
      I1 = index1*3; 
      I2 = index2*3; 
      EnergyNonbond enb;
/// Atom-vector in AtomTable is empty
//      enb._Atom1 = this->_AtomTable->elt_atom(index1);
//      enb._Atom2 = this->_AtomTable->elt_atom(index2);
      enb.term.dA = dA;
      enb.term.dC = dC;
      enb._Charge1 = charge11;
      enb._Charge2 = charge22;
      enb.term.I1 = I1;
      enb.term.I2 = I2;
      if (termIndex<count) {
        this->_Terms[termIndex] = enb;
      } else {
        SIMPLE_ERROR(BF("Overflowed the _Terms array with termIndex=%lu and count = %lu\n") % termIndex % count);
      }
      ++termIndex;
//      printf( "nonbond index1 name %s index2 %s\n",  this->_AtomTable->elt_atom_name(index1), this->_AtomTable->elt_atom_name(index2));
#ifdef DEBUG_NONBOND_TERM
      if ( this->_DebugEnergy ) {
        std::string key;
        std::string atom1Name = gc::As<core::Symbol_sp>((*this->_atom_name_vector)[index1])->symbolNameAsString();
        std::string atom2Name = gc::As<core::Symbol_sp>((*this->_atom_name_vector)[index2])->symbolNameAsString();
        if ( atom1Name < atom2Name ) {
          key = atom1Name+"-"+atom2Name;
        } else {
          key = atom2Name+"-"+atom1Name;
        }
        //printf( "nonbond %s  vdw(%lf) electrostatic(%lf)\n", (double)this->_EnergyVdw,  this->_EnergyElectrostatic );
        LOG_ENERGY(BF( "MEISTER nonbond %s args cando\n")% key );
        LOG_ENERGY(BF( "MEISTER nonbond %s dA %5.3lf\n")% key % dA );
        LOG_ENERGY(BF( "MEISTER nonbond %s dC %5.3lf\n")% key % dC );
        LOG_ENERGY(BF( "MEISTER nonbond %s dQ1Q2 %5.3lf\n")% key % dQ1Q2 );
        LOG_ENERGY(BF( "MEISTER nonbond %s x1 %5.3lf %d\n")% key % x1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s y1 %5.3lf %d\n")% key % y1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s z1 %5.3lf %d\n")% key % z1 % (I1/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s x2 %5.3lf %d\n")% key % x2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s y2 %5.3lf %d\n")% key % y2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s z2 %5.3lf %d\n")% key % z2 % (I2/3+1) );
        LOG_ENERGY(BF( "MEISTER nonbond %s results\n")% key );
        LOG_ENERGY(BF( "MEISTER nonbond %s evdw %lf\n")% key % Evdw);
        LOG_ENERGY(BF( "MEISTER nonbond %s eeel %lf\n")% key % Eeel);
        LOG_ENERGY(BF( "MEISTER nonbond %s Enonbond(evdw+eeel) %lf\n")% key % (Evdw+Eeel) );
//        if ( calcForce ) { 
//          LOG_ENERGY(BF( "MEISTER nonbond %s fx1 %lf %d\n")% key % fx1 % (I1/3+1) );
//          LOG_ENERGY(BF( "MEISTER nonbond %s fy1 %lf %d\n")% key % fy1 % (I1/3+1) );
//          LOG_ENERGY(BF( "MEISTER nonbond %s fz1 %lf %d\n")% key % fz1 % (I1/3+1) );
//          LOG_ENERGY(BF( "MEISTER nonbond %s fx2 %lf %d\n")% key % fx2 % (I2/3+1) );
//          LOG_ENERGY(BF( "MEISTER nonbond %s fy2 %lf %d\n")% key % fy2 % (I2/3+1) );
//          LOG_ENERGY(BF( "MEISTER nonbond %s fz2 %lf %d\n")% key % fz2 % (I2/3+1) );
//        }
        LOG_ENERGY(BF( "MEISTER nonbond %s stop\n")% key );
      }
#endif
    }
  }
  printf("%s:%d:%s    values -> %lu\n", __FILE__, __LINE__, __FUNCTION__, this->_Terms.size());

//  printf( "Nonbond energy vdw(%lf) electrostatic(%lf)\n", (double)this->_EnergyVdw,  this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy }\n"));
}










void	EnergyNonbond_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
                                                                                NVector_sp 	pos)
{_OF();
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
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
    double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
    int	I1, I2,i;
    gctools::Vec0<EnergyNonbond>::iterator nbi;
    for ( i=0,nbi=this->_Terms.begin();
          nbi!=this->_Terms.end(); nbi++,i++ ) {
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
      int index = i;
#include <cando/chem/energy_functions/_Nonbond_debugFiniteDifference.cc>

    }
  }
}



int	EnergyNonbond_O::checkForBeyondThresholdInteractions(
                                                             stringstream& info, NVector_sp pos )
{_OF();
  int	fails = 0;

  this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------
#undef NONBOND_CALC_FORCE
#undef NONBOND_CALC_DIAGONAL_HESSIAN
#undef NONBOND_CALC_OFF_DIAGONAL_HESSIAN
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

  if ( this->isEnabled() )
  {
    _BLOCK_TRACE("NonbondEnergy finiteDifference comparison");
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
    double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
    int	I1, I2,i;
    gctools::Vec0<EnergyNonbond>::iterator nbi;
    for ( i=0,nbi=this->_Terms.begin();
          nbi!=this->_Terms.end(); nbi++,i++ ) 
    {
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
      if ( NonbondDistance < this->_ErrorThreshold ) {
        Atom_sp a1, a2;
        a1 = (*nbi)._Atom1;
        a2 = (*nbi)._Atom2;
        info<< "NonbondDeviation";
//		info<< a1->getAbsoluteIdPath() << " ";
//		info<< a2->getAbsoluteIdPath() << " ";
        info<< "value " << NonbondDistance << " ";
        info<<"threshold " << this->_ErrorThreshold;
        info << a1->getName() << " ";
        info << a2->getName() << " ";
        info << std::endl;
        this->_BeyondThresholdTerms.push_back(*nbi);
        fails++;
      }
    }
  }
  return fails;
}

void EnergyNonbond_O::initialize()
{
  this->Base::initialize();
  this->setErrorThreshold(1.0);
  this->setDielectricConstant(80.0);
  this->setVdwScale(1.0);
  this->setElectrostaticScale(1.0);
}



void EnergyNonbond_O::addTerm(const EnergyNonbond& term)
{
  this->_Terms.push_back(term);
}

void EnergyNonbond_O::fields(core::Record_sp node)
{
node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

string EnergyNonbond_O::beyondThresholdInteractionsAsString()
{
  return component_beyondThresholdInteractionsAsString<EnergyNonbond_O,EnergyNonbond>(*this);
}


void EnergyNonbond_O::construct14InteractionTerms(AtomTable_sp atomTable, Matter_sp matter, core::T_sp forceField, core::T_sp activeAtoms, bool show_progress)
{
  _OF();
  {
    EnergyNonbond energyNonbond;
    Loop loop;
    loop.loopTopGoal(matter,PROPERS);
    size_t terms(0);
    while ( loop.advanceLoopAndProcess() ) {
      Atom_sp a1 = loop.getAtom1();
      Atom_sp a4 = loop.getAtom4();
      if ( activeAtoms.notnilp() &&
           ( !inAtomSet(activeAtoms,a1)
             || !inAtomSet(activeAtoms,a4) )) continue;
      auto ea1 = atomTable->getEnergyAtomPointer(a1);
      auto ea4 = atomTable->getEnergyAtomPointer(a4);
      energyNonbond.defineFrom(forceField,true,&*ea1,&*ea4,this->asSmartPtr());
      LOG(BF("About to addTerm"));
      this->addTerm(energyNonbond);
      LOG(BF("Returned from addTerm"));
      ++terms;
    }
    if (show_progress) core::write_bf_stream(BF("Built 14 interaction table with %d terms\n") % terms);
  }
}

void EnergyNonbond_O::constructNonbondTermsFromAtomTable(bool ignore14s, AtomTable_sp atomTable, core::T_sp nbForceField, bool show_progress)
{
      printf("%s:%d In :constructNonbondTermsFromAtomTable\n", __FILE__, __LINE__ );

  // ------------------------------------------------------------
  //
  // The old code created terms for each nonbonded interaction
  // this will use waaaay too much memory for large systems
  //
  this->_UsesExcludedAtoms = false;
  if ( atomTable->getNumberOfAtoms() > 2 ) {
    LOG(BF("Defining NONBONDS") );
    gctools::Vec0<EnergyAtom>::iterator	iea1;
    gctools::Vec0<EnergyAtom>::iterator	iea2;
    size_t total_comparisons = atomTable->getNumberOfAtoms()*atomTable->getNumberOfAtoms()/2;
    if (show_progress) {
      core::write_bf_stream(BF("For nonbonded interactions, about to carry out %zu atom-to-atom comparisons\n") % total_comparisons);
    }
    for ( iea1 = atomTable->begin();
          iea1 != atomTable->end()-1; iea1++ ) { //Was iea1 != atomTable->end()-1
      for ( iea2 = iea1+1; iea2 != atomTable->end(); iea2++ ) { // Was iea2 != atomTable->end()
        if (!iea1->inBondOrAngle(iea2->atom())) {
          bool in14 = iea1->relatedBy14(iea2->atom());
          // 14s are added separately to the Terms
          if (in14 ) {
            if ( !ignore14s ) {
 //             LOG_ENERGY(BF("Nonbonded interaction between %s - %s in14[%d]\n") % _rep_(iea1->atom()) % _rep_(iea2->atom()) % in14 );
              EnergyNonbond energyNonbond;
              if ( energyNonbond.defineFrom(nbForceField, in14,
                                            &(*iea1),&(*iea2),this->sharedThis<EnergyNonbond_O>()) )  {
                this->addTerm(energyNonbond);
                LOG_ENERGY(BF("nonbond  interaction between %s - %s in14[%d] dA %lf\n") % _rep_(iea1->atom()) % _rep_(iea2->atom()) % in14  % energyNonbond.term.dA);                 
              }
            }
          } else {
//            LOG_ENERGY(BF("Nonbonded interaction between %s - %s in14[%d]\n") % _rep_(iea1->atom()) % _rep_(iea2->atom()) % in14 );
            EnergyNonbond energyNonbond;
            if ( energyNonbond.defineFrom(nbForceField, in14,
                                          &(*iea1),&(*iea2),this->sharedThis<EnergyNonbond_O>()) )  {
              this->addTerm(energyNonbond);
              LOG_ENERGY(BF("nonbond  interaction between %s - %s in14[%d] dA %lf\n") % _rep_(iea1->atom()) % _rep_(iea2->atom()) % in14  % energyNonbond.term.dA);                 

            }
          }
        }
      }
    }
  } else {
    LOG_ENERGY(BF("There are no non-bonds\n"));
  }
}

void EnergyNonbond_O::constructExcludedAtomListFromAtomTable(AtomTable_sp atomTable, core::T_sp nbForceField, bool show_progress)
{
//  printf("%s:%d In :constructExcludedATomListFromAtomTable\n", __FILE__, __LINE__ );
    // ------------------------------------------------------------
  //
  // The new code creates an excluded atom list
  //
  this->_UsesExcludedAtoms = true;
  this->_AtomTable = atomTable;
  this->_FFNonbondDb = nbForceField;
  core::T_mv values_mv = this->_AtomTable->calculate_excluded_atom_list();
  core::SimpleVector_int32_t_sp number_of_excluded_atoms = gc::As<core::SimpleVector_int32_t_sp>(values_mv);
  core::SimpleVector_int32_t_sp excluded_atoms_list = gc::As<core::SimpleVector_int32_t_sp>(values_mv.second());
  this->_NumberOfExcludedAtomIndices = number_of_excluded_atoms;
  this->_ExcludedAtomIndices = excluded_atoms_list;
}

SYMBOL_EXPORT_SC_(KeywordPkg,da);
SYMBOL_EXPORT_SC_(KeywordPkg,dc);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,charge1);
SYMBOL_EXPORT_SC_(KeywordPkg,charge2);
SYMBOL_EXPORT_SC_(KeywordPkg,atomname1);
SYMBOL_EXPORT_SC_(KeywordPkg,atomname2);

CL_DEFMETHOD core::List_sp EnergyNonbond_O::extract_vectors_as_alist() const{
  size_t size = this->_Terms.size();
  core::SimpleVector_double_sp da_vec = core::SimpleVector_double_O::make(size);
  core::SimpleVector_double_sp dc_vec = core::SimpleVector_double_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_double_sp charge1_vec = core::SimpleVector_double_O::make(size);
  core::SimpleVector_double_sp charge2_vec = core::SimpleVector_double_O::make(size);
  core::SimpleVector_sp atomname1_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atomname2_vec = core::SimpleVector_O::make(size);
  for (size_t i=0; i<size; ++i){
    const EnergyNonbond& entry = this->_Terms[i];
    (*da_vec)[i] = entry.term.dA;
    (*dc_vec)[i] = entry.term.dC;
    (*i1_vec)[i] = entry.term.I1;
    (*i2_vec)[i] = entry.term.I2;
    (*charge1_vec)[i] = this->_AtomTable->elt_charge(entry.term.I1/3);
    (*charge2_vec)[i] = this->_AtomTable->elt_charge(entry.term.I2/3);
    (*atomname1_vec)[i] = this->_AtomTable->elt_atom_name(entry.term.I1/3);
    (*atomname2_vec)[i] = this->_AtomTable->elt_atom_name(entry.term.I2/3);
  }
  return core::Cons_O::createList(core::Cons_O::create(kw::_sym_da, da_vec),
                                  core::Cons_O::create(kw::_sym_dc, dc_vec),
                                  core::Cons_O::create(kw::_sym_i1, i1_vec),
                                  core::Cons_O::create(kw::_sym_i2, i2_vec),
                                  core::Cons_O::create(kw::_sym_charge1, charge1_vec),
                                  core::Cons_O::create(kw::_sym_charge2, charge2_vec),
                                  core::Cons_O::create(kw::_sym_atomname1 ,atomname1_vec),
                                  core::Cons_O::create(kw::_sym_atomname2 ,atomname2_vec));
}


CL_DEFMETHOD core::T_sp EnergyNonbond_O::getFFNonbondDb() {
  if (this->_FFNonbondDb.boundp()) return this->_FFNonbondDb;
  SIMPLE_ERROR(BF("The _FFNonbondDb of an EnergyNonbond has not been initialized"));
}



SYMBOL_EXPORT_SC_(KeywordPkg,ntypes);
SYMBOL_EXPORT_SC_(KeywordPkg,atom_name_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,charge_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,mass_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,atomic_number_vector);
SYMBOL_EXPORT_SC_(KeywordPkg,ico_vec);
SYMBOL_EXPORT_SC_(KeywordPkg,iac_vec);
SYMBOL_EXPORT_SC_(KeywordPkg,local_typej_vec);
SYMBOL_EXPORT_SC_(KeywordPkg,cn1_vec);
SYMBOL_EXPORT_SC_(KeywordPkg,cn2_vec);
SYMBOL_EXPORT_SC_(KeywordPkg,excluded_atoms_list);
SYMBOL_EXPORT_SC_(KeywordPkg,number_excluded_atoms);
SYMBOL_EXPORT_SC_(KeywordPkg,atom_table);



CL_DEFMETHOD void EnergyNonbond_O::constructNonbondTermsFromAList(core::List_sp values)
{
//  printf("%s:%d:%s    values -> %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(values).c_str());
  this->_ntypes =               translate::from_object<size_t>(safe_alist_lookup<core::T_sp>(values,kw::_sym_ntypes))._v;          // ntypes
  this->_atom_name_vector =     safe_alist_lookup<core::SimpleVector_sp>(values,kw::_sym_atom_name_vector);  // atom-name-vector
  this->_charge_vector =        safe_alist_lookup<core::SimpleVector_double_sp>(values,kw::_sym_charge_vector);          // charge-vector
  this->_mass_vector =          safe_alist_lookup<core::SimpleVector_double_sp>(values,kw::_sym_mass_vector);            // masses
  this->_atomic_number_vector = safe_alist_lookup<core::SimpleVector_int32_t_sp>(values,kw::_sym_atomic_number_vector);    // vec
  this->_ico_vec =              safe_alist_lookup<core::SimpleVector_int32_t_sp>(values,kw::_sym_ico_vec);             // ico-vec
  this->_iac_vec =              safe_alist_lookup<core::SimpleVector_int32_t_sp>(values,kw::_sym_iac_vec);             // iac-vec
//  this->_local_typej_vec =      safe_alist_lookup(values,kw::_sym_local_typej_vec);      // local-typej-vec
  this->_cn1_vec =              safe_alist_lookup<core::SimpleVector_double_sp>(values,kw::_sym_cn1_vec);
  this->_cn2_vec =              safe_alist_lookup<core::SimpleVector_double_sp>(values,kw::_sym_cn2_vec);
}

CL_DEFMETHOD void EnergyNonbond_O::setNonbondExcludedAtomInfo(AtomTable_sp atom_table, core::SimpleVector_int32_t_sp excluded_atoms_list, core::SimpleVector_int32_t_sp number_excluded_atoms) {
  this->_AtomTable = atom_table;
  this->_ExcludedAtomIndices = excluded_atoms_list;
  this->_NumberOfExcludedAtomIndices = number_excluded_atoms;
    printf("%s:%d:%s   Exiting\n", __FILE__, __LINE__, __FUNCTION__ );

}

//core::List_sp termAtIndex(size_t index) const
//{


//}

};

