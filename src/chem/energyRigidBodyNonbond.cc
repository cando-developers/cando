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
#include <cando/chem/energyRigidBodyNonbond.h>
#include <cando/chem/rigidBodyEnergyFunction.h>
#include <clasp/core/profiler.h>
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
#include <cando/chem/nVector.h>
#include <cando/chem/loop.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>


#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
//#define LOG_ENERGY BFORMAT_T

namespace chem
{

#ifdef XML_ARCHIVE
void	EnergyRigidBodyNonbond::archive(core::ArchiveP node)
{
  node->attribute("_Is14",this->_Is14);
  node->attribute("_A",this->_A);
  node->attribute("_C",this->_C);
  node->attribute("_Charge1",this->_Charge1);
  node->attribute("_Charge2",this->_Charge2);
  node->attribute("_RStar",this->_RStar);
  node->attribute("dQ1Q2",this->term.dQ1Q2);
  node->attribute("dA",this->term.dA);
  node->attribute("dC",this->term.dC);
  node->attribute("I1",this->term.I1);
  node->attribute("I2",this->term.I2);
  node->attribute("a1",this->_Atom1);
  node->attribute("a2",this->_Atom2);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
  node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
  node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
  node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include <cando/energy-functions/_NONBONDRB_debugEvalSerialize.cc>
#endif //]
}
#endif

#if 0
//
// Return true if we could fill the energyNonbond term
// otherwise false usually if we don't recognize one of the atom types like DU
//
bool	EnergyRigidBodyNonbond::defineFrom(ForceField_sp	forceField,
                                  bool		is14,
                                  EnergyAtom	*iea1,
                                  EnergyAtom	*iea2,
                                  EnergyRigidBodyNonbond_sp energyNonbond)
{_OF();
  LOG(BF("defineFrom"));
  FFNonbond_sp				ffNonbond1;
  FFNonbond_sp				ffNonbond2;
  double				epsilonij;
  double				vdwScale;
  double				electrostaticScale;
  this->_Is14 = is14;
  this->_Atom1 = iea1->atom();
  this->_Atom2 = iea2->atom();
  core::Symbol_sp t1 = iea1->atom()->getType();
  core::Symbol_sp t2 = iea2->atom()->getType();
  LOG(BF("Defining nonbond between types: %s - %s") % _rep_(t1) % _rep_(t2));
  ASSERT(forceField->_Nonbonds&&forceField->_Nonbonds.notnilp());
  LOG(BF("forceField->_Nonbonds @%p   .notnilp()->%d") % forceField->_Nonbonds.raw_() % forceField->_Nonbonds.notnilp());
  ffNonbond1 = forceField->_Nonbonds->findType(t1);
  ffNonbond2 = forceField->_Nonbonds->findType(t2);
  ANN(ffNonbond1);
  if ( ffNonbond1.nilp() )
  {
//     	SIMPLE_ERROR(BF("Unknown force field type(",iea1->_Atom->getType().c_str(),") for non-bonded interaction"));
    return false;
  }
  ANN(ffNonbond2);
  if ( ffNonbond2.nilp() )
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
#endif

#if 0
double	EnergyRigidBodyNonbond::getDistance()
{
  Vector3	pos1, pos2;
  pos1 = this->_Atom1->getPosition();
  pos2 = this->_Atom2->getPosition();
  return geom::calculateDistance(pos1,pos2);
}
#endif

#if 0
adapt::QDomNode_sp	EnergyRigidBodyNonbond::asXml()
{
  adapt::QDomNode_sp	node;
  Vector3	vdiff;

  node = adapt::QDomNode_O::create(lisp,"EnergyRigidBodyNonbond");
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
#include <_NONBONDRB_debugEvalXml.cc>
  node->addChild(xml);
#endif
  node->addAttributeDoubleScientific("dA",this->term.dA);
  node->addAttributeDoubleScientific("dC",this->term.dC);
  node->addAttributeDoubleScientific("dQ1Q2",this->term.dQ1Q2);
  return node;
}

void	EnergyRigidBodyNonbond::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
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



#if 0
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
#include <cando/energy-functions/_NONBONDRB_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/energy-functions/_NONBONDRB_termCode.cc>

  return Energy;
}
#endif


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


#if 0
void	EnergyRigidBodyNonbond_O::dumpTerms()
{
  gctools::Vec0<EnergyRigidBodyNonbond>::iterator	eni;
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
#endif


void	EnergyRigidBodyNonbond_O::setupHessianPreconditioner(
                                                    NVector_sp nvPosition,
                                                    AbstractLargeSquareMatrix_sp m )
{
  SIMPLE_ERROR(BF("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!!"));
}

double	EnergyRigidBodyNonbond_O::evaluateAll(NVector_sp 	pos,
                                     bool 		calcForce,
                                     gc::Nilable<NVector_sp> 	force,
                                     bool		calcDiagonalHessian,
                                     bool		calcOffDiagonalHessian,
                                     gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                     gc::Nilable<NVector_sp>	hdvec, 
                                     gc::Nilable<NVector_sp> 	dvec )
{
  double vdwScale = this->getVdwScale();
  double electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
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
#define	NONBONDRB_EEEL_ENERGY_ACCUMULATE(e) {this->_EnergyElectrostatic +=(e);}
#undef	NONBONDRB_EVDW_ENERGY_ACCUMULATE
#define	NONBONDRB_EVDW_ENERGY_ACCUMULATE(e) {this->_EnergyVdw+=(e);}
#undef	NONBONDRB_ENERGY_ACCUMULATE
#define	NONBONDRB_ENERGY_ACCUMULATE(e) {};
#undef	NONBONDRB_FORCE_ACCUMULATE
#undef	NONBONDRB_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBONDRB_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBONDRB_FORCE_ACCUMULATE 		ForceAcc
#define	NONBONDRB_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBONDRB_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
  if ( !this->isEnabled() ) return 0.0;
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
  LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_NONBONDRB_termDeclares.cc>
#pragma clang diagnostic pop
  double dA,dC,dQ1Q2;
  double xh1, yh1, zh1;
  double xh2, yh2, zh2;
  double ak, bk, ck, dk, xk, yk, zk;
  double al, bl, cl, dl, xl, yl, zl;
  int	I1, I2;
  size_t I1start = 0;
  for ( size_t iI1 = 0; iI1<(this->_RigidBodyEndAtom->length()-1); ++iI1 ) {
    size_t I1end = (*this->_RigidBodyEndAtom)[iI1];
    size_t I2start = I1end;
    for (size_t iI2 = iI1+1; iI2<this->_RigidBodyEndAtom->length(); ++iI2 ) {
      size_t I2end = (*this->_RigidBodyEndAtom)[iI2];
      for ( size_t I1cur = I1start; I1cur<I1end; ++I1cur ) {
        RigidBodyAtomInfo& ea1 = this->_AtomInfoTable[I1cur];
        double charge1 = ea1._Charge;
        double electrostatic_scaled_charge1 = charge1*electrostaticScale;
        for ( size_t I2cur = I2start; I2cur<I2end; ++I2cur ) {
          RigidBodyAtomInfo& ea2 = this->_AtomInfoTable[I2cur];
          double rStar = ea1._Radius+ea2._Radius;
          double epsilonij = sqrt(ea1._Epsilon*ea2._Epsilon);
          double rStar2 = rStar*rStar;
          double rStar6 = rStar2*rStar2*rStar2;
          double rStar12 = rStar6*rStar6;
          dA = epsilonij*rStar12*vdwScale;
          dC = 2.0*epsilonij*rStar6*vdwScale;
          double charge2 = ea2._Charge;
          dQ1Q2 = electrostatic_scaled_charge1*charge2;
          I1 = iI1*7;
          I2 = iI2*7;
#include <cando/energy-functions/_NONBONDRB_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
          nbi->_calcForce = calcForce;
          nbi->_calcDiagonalHessian = calcDiagonalHessian;
          nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/energy-functions/_NONBONDRB_debugEvalSet.cc>
#endif //]
#ifdef DEBUG_NONBOND_TERM
          if ( this->_DebugEnergy ) {
            std::string key;
            if ( ea1._Atom->getName()->symbolNameAsString()
                 < ea2._Atom->getName()->symbolNameAsString() ) {
              key = ea1._Atom->getName()->symbolNameAsString()+"-"+ea2._Atom->getName()->symbolNameAsString();
            } else {
              key = ea2._Atom->getName()->symbolNameAsString()+"-"+ea1._Atom->getName()->symbolNameAsString();
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
    }
  }
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF( "Nonbond energy }\n"));
  return this->_TotalEnergy;
}
    
    



#if 0
ALWAYS_INLINE void _calculate_nonbond_term(EnergyRigidBodyNonbond_O& me, double vdwScale, double electrostatic_scaled_charge1, FFNonbond_sp ffNonbond1, FFNonbond_sp ffNonbond2, AtomTable_O::iterator iea1, AtomTable_O::iterator iea2,
                                           NVector_sp 	pos,
                                           bool 		calcForce,
                                           gc::Nilable<NVector_sp> 	force,
                                           bool		calcDiagonalHessian,
                                           bool		calcOffDiagonalHessian,
                                           gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                           gc::Nilable<NVector_sp>	hdvec, 
                                           gc::Nilable<NVector_sp> 	dvec )
{
#define NONBOND_CALC_FORCE
#define NONBOND_CALC_DIAGONAL_HESSIAN
#define NONBOND_CALC_OFF_DIAGONAL_HESSIAN
#undef	NONBOND_SET_PARAMETER
#define	NONBOND_SET_PARAMETER(x)	{}
#undef	NONBOND_SET_POSITION
#define	NONBOND_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef	NONBOND_EEEL_ENERGY_ACCUMULATE
#define	NONBOND_EEEL_ENERGY_ACCUMULATE(e) {me._EnergyElectrostatic +=(e);}
#undef	NONBOND_EVDW_ENERGY_ACCUMULATE
#define	NONBOND_EVDW_ENERGY_ACCUMULATE(e) {me._EnergyVdw+=(e);}
#undef	NONBOND_ENERGY_ACCUMULATE
#define	NONBOND_ENERGY_ACCUMULATE(e) {};
#undef	NONBOND_FORCE_ACCUMULATE
#undef	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#undef	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	NONBOND_FORCE_ACCUMULATE 		ForceAcc
#define	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/energy-functions/_NONBONDRB_termDeclares.cc>
#pragma clang diagnostic pop
  double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
  double rStar = ffNonbond1->getRadius_Angstroms()+ffNonbond2->getRadius_Angstroms();
  // OPTIMIZE: Can I precalculate the square-roots!!!!!!!!
  double epsilonij = sqrt(ffNonbond1->getEpsilon_kCal()*ffNonbond2->getEpsilon_kCal()); 
  double rStar6 = pow(rStar,6.0);
  double rStar12 = rStar*rStar;
  dA = epsilonij*rStar12*vdwScale;
  dC = 2.0*epsilonij*rStar6*vdwScale;
  double charge2 = iea2->atom()->getCharge();
  dQ1Q2 = electrostatic_scaled_charge1*charge2;
  LOG(BF( "Calc dQ1Q2 electrostaticScale= %lf")% (double)(electrostaticScale));
  LOG(BF( "Calc dQ1Q2 Dielectric constant = %lf")% (double)(dielectricConstant));
  LOG(BF( "Calc dQ1Q2 Charge1 = %lf")% (double)(me._Charge1));
  LOG(BF( "Calc dQ1Q2 Charge2 = %lf")% (double)(me._Charge2));
  LOG(BF( "dQ1Q2 = %lf")% (double)(me.term.dQ1Q2));
  int I1 = iea1->coordinateIndexTimes3();
  int I2 = iea2->coordinateIndexTimes3();
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
  if ( me._Debug_NumberOfNonbondTermsToCalculate > 0 ) {
    if ( i>= me._Debug_NumberOfNonbondTermsToCalculate ) {
      break;
    }
  }
#endif
#include <cando/energy-functions/_NONBONDRB_termCode.cc>
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
  nbi->_calcForce = calcForce;
  nbi->_calcDiagonalHessian = calcDiagonalHessian;
  nbi->_calcOffDiagonalHessian = calcOffDiagonalHessian;
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ nbi->eval.var=val;};
#include <cando/energy-functions/_NONBONDRB_debugEvalSet.cc>
#endif //]
#ifdef DEBUG_NONBOND_TERM
  if ( me._DebugEnergy ) {
    LOG_ENERGY(BF( "MEISTER nonbond %d args cando\n")% (i+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d dA %5.3lf\n")% (i+1) % (nbi->term.dA) );
    LOG_ENERGY(BF( "MEISTER nonbond %d dC %5.3lf\n")% (i+1) % (nbi->term.dC) );
    LOG_ENERGY(BF( "MEISTER nonbond %d dQ1Q2 %5.3lf\n")% (i+1) % (nbi->term.dQ1Q2) );
    LOG_ENERGY(BF( "MEISTER nonbond %d x1 %5.3lf %d\n")% (i+1) % x1 % (I1/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d y1 %5.3lf %d\n")% (i+1) % y1 % (I1/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d z1 %5.3lf %d\n")% (i+1) % z1 % (I1/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d x2 %5.3lf %d\n")% (i+1) % x2 % (I2/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d y2 %5.3lf %d\n")% (i+1) % y2 % (I2/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d z2 %5.3lf %d\n")% (i+1) % z2 % (I2/3+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d results\n")% (i+1) );
    LOG_ENERGY(BF( "MEISTER nonbond %d Enonbond %lf\n")% (i+1) % (Evdw+Eeel) );
    LOG_ENERGY(BF( "MEISTER nonbond %d evdw %lf\n")% (i+1) % Evdw);
    LOG_ENERGY(BF( "MEISTER nonbond %d eeel %lf\n")% (i+1) % Eeel);
    if ( calcForce ) {
      LOG_ENERGY(BF( "MEISTER nonbond %d fx1 %lf %d\n")% (i+1) % fx1 % (I1/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fy1 %lf %d\n")% (i+1) % fy1 % (I1/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fz1 %lf %d\n")% (i+1) % fz1 % (I1/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fx2 %lf %d\n")% (i+1) % fx2 % (I2/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fy2 %lf %d\n")% (i+1) % fy2 % (I2/3+1) );
      LOG_ENERGY(BF( "MEISTER nonbond %d fz2 %lf %d\n")% (i+1) % fz2 % (I2/3+1) );
    }
    LOG(BF( "MEISTER nonbond %d stop\n")% (i+1) );
  }
#endif
}
#endif






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
#include <cando/energy-functions/_NONBONDRB_termDeclares.cc>
#pragma clang diagnostic pop
    double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
    int	I1, I2,i;
    gctools::Vec0<EnergyRigidBodyNonbond>::iterator nbi;
    for ( i=0,nbi=this->_Terms.begin();
          nbi!=this->_Terms.end(); nbi++,i++ ) {
#include <cando/energy-functions/_NONBONDRB_termCode.cc>
      int index = i;
#include <cando/energy-functions/_NONBONDRB_debugFiniteDifference.cc>

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
