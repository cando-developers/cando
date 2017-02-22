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
#define	DEBUG_LEVEL_NONE

#include <clasp/core/foundation.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
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
//#define LOG_ENERGY BFORMAT_T

namespace chem
{

#ifdef XML_ARCHIVE
void	EnergyNonbond::archive(core::ArchiveP node)
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
#include <cando/chem/energy_functions/_Nonbond_debugEvalSerialize.cc>
#endif //]
}
#endif

//
// Return true if we could fill the energyNonbond term
// otherwise false usually if we don't recognize one of the atom types like DU
//
bool	EnergyNonbond::defineFrom(ForceField_sp	forceField,
                                  bool		is14,
                                  EnergyAtom	*iea1,
                                  EnergyAtom	*iea2,
                                  EnergyNonbond_sp energyNonbond)
{
  FFNonbond_sp				ffNonbond1;
  FFNonbond_sp				ffNonbond2;
  double				epsilonij;
  double				vdwScale;
  double				electrostaticScale;
  this->_Is14 = is14;
  this->_Atom1 = iea1->atom();
  this->_Atom2 = iea2->atom();
  ffNonbond1 = forceField->_Nonbonds->findType(iea1->atom()->getType());
  ffNonbond2 = forceField->_Nonbonds->findType(iea2->atom()->getType());
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

void	EnergyNonbond_O::evaluateAll(NVector_sp 	pos,
                                     bool 		calcForce,
                                     gc::Nilable<NVector_sp> 	force,
                                     bool		calcDiagonalHessian,
                                     bool		calcOffDiagonalHessian,
                                     gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                     gc::Nilable<NVector_sp>	hdvec, 
                                     gc::Nilable<NVector_sp> 	dvec )
{
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
}
    
    

void	EnergyNonbond_O::evaluateTerms(NVector_sp 	pos,
                                     bool 		calcForce,
                                     gc::Nilable<NVector_sp> 	force,
                                     bool		calcDiagonalHessian,
                                     bool		calcOffDiagonalHessian,
                                     gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                     gc::Nilable<NVector_sp>	hdvec, 
                                     gc::Nilable<NVector_sp> 	dvec )
{
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
  LOG(BF( "Nonbond energy vdw(%lf) electrostatic(%lf)\n")% (double)this->_EnergyVdw % this->_EnergyElectrostatic );
  LOG(BF("Nonbond energy }\n"));
}



#if 0
ALWAYS_INLINE void _calculate_nonbond_term(EnergyNonbond_O& me, double vdwScale, double electrostatic_scaled_charge1, FFNonbond_sp ffNonbond1, FFNonbond_sp ffNonbond2, AtomTable_O::iterator iea1, AtomTable_O::iterator iea2,
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
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
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
                                          


void	EnergyNonbond_O::evaluateUsingExcludedAtoms(NVector_sp 	pos,
                                                    bool 		calcForce,
                                                    gc::Nilable<NVector_sp> 	force,
                                                    bool		calcDiagonalHessian,
                                                    bool		calcOffDiagonalHessian,
                                                    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                    gc::Nilable<NVector_sp>	hdvec, 
                                                    gc::Nilable<NVector_sp> 	dvec )
{
  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  core::NativeVector_int_sp numberOfExcludedAtoms = this->_NumberOfExcludedAtomIndices;
  core::NativeVector_int_sp excludedAtomIndices = this->_ExcludedAtomIndices;
  double vdwScale = this->getVdwScale();
  double electrostaticScale = this->getElectrostaticScale()*ELECTROSTATIC_MODIFIER/this->getDielectricConstant();
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
  if ( !this->isEnabled() ) return;
	    // If you are going to use openmp here, you need to control access to the force and hessian
	    // arrays so that only one thread updates each element at a time.
  LOG(BF("Nonbond component is enabled") );
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
  double x1,y1,z1,x2,y2,z2,dA,dC,dQ1Q2;
  int	I1, I2;
  int i = 0;
  int maxIndex = this->_AtomTable->getNumberOfAtoms()*3;
  auto iea1_end = this->_AtomTable->end()-1;
  int index1 = 0;
  int excludedAtomIndex = 0;
  for ( auto iea1 = this->_AtomTable->begin(); iea1 != iea1_end; ++iea1, ++index1 ) {
    LOG(BF("%s ====== top of outer loop - index1 = %d\n") % __FUNCTION__ % index1 );
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    LOG(BF("Read numberOfExcludedAtomsRemaining from numberOfExcludedAtoms[%d]= %d\n") % index1 % numberOfExcludedAtomsRemaining);
          // Skip 0 in excluded atom list that amber requires
    if (numberOfExcludedAtomsRemaining<0) {
      ++excludedAtomIndex;
    }
    int index2 = index1+1;
    FFNonbond_sp ffNonbond1 = this->_FFNonbondDb->getFFNonbondUsingTypeIndex(iea1->_TypeIndex);
    double charge1 = iea1->atom()->getCharge();
    double electrostatic_scaled_charge1 = charge1*electrostaticScale;
    for ( auto iea2 = iea1+1; iea2 != this->_AtomTable->end(); iea2++, ++i, ++index2 ) {
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
      FFNonbond_sp ffNonbond2 = this->_FFNonbondDb->getFFNonbondUsingTypeIndex(iea2->_TypeIndex);
      double rStar = ffNonbond1->getRadius_Angstroms()+ffNonbond2->getRadius_Angstroms();
      double epsilonij = sqrt(ffNonbond1->getEpsilon_kCal()*ffNonbond2->getEpsilon_kCal());
      double rStar2 = rStar*rStar;
      double rStar6 = rStar2*rStar2*rStar2;
      double rStar12 = rStar6*rStar6;
      dA = epsilonij*rStar12*vdwScale;
      dC = 2.0*epsilonij*rStar6*vdwScale;
      double charge2 = iea2->atom()->getCharge();
      dQ1Q2 = electrostatic_scaled_charge1*charge2;
      I1 = iea1->coordinateIndexTimes3();
      I2 = iea2->coordinateIndexTimes3();
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
        if ( iea1->_AtomName->symbolNameAsString()
             < iea2->_AtomName->symbolNameAsString() ) {
          key = iea1->_AtomName->symbolNameAsString()+"-"+iea2->_AtomName->symbolNameAsString();
        } else {
          key = iea2->_AtomName->symbolNameAsString()+"-"+iea1->_AtomName->symbolNameAsString();
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


#ifdef XML_ARCHIVE
void EnergyNonbond_O::archiveBase(core::ArchiveP node)
{
  this->Base::archiveBase(node);
  archiveEnergyComponentTerms<EnergyNonbond_O,EnergyNonbond>(node,*this);
}
#endif

string EnergyNonbond_O::beyondThresholdInteractionsAsString()
{
  return component_beyondThresholdInteractionsAsString<EnergyNonbond_O,EnergyNonbond>(*this);
}


void EnergyNonbond_O::construct14InteractionTerms(AtomTable_sp atomTable, Matter_sp matter, ForceField_sp forceField, core::T_sp activeAtoms, bool show_progress)
{
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
      this->addTerm(energyNonbond);
      ++terms;
    }
    if (show_progress) BFORMAT_T(BF("Built 14 interaction table with %d terms\n") % terms);
  }
}

void EnergyNonbond_O::constructNonbondTermsFromAtomTable(bool ignore14s, AtomTable_sp atomTable, ForceField_sp forceField, bool show_progress)
{
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
      BFORMAT_T(BF("For nonbonded interactions, about to carry out %zu atom-to-atom comparisons\n") % total_comparisons);
    }
    for ( iea1 = atomTable->begin();
          iea1 != atomTable->end()-1; iea1++ ) { //Was iea1 != atomTable->end()-1
      for ( iea2 = iea1+1; iea2 != atomTable->end(); iea2++ ) { // Was iea2 != atomTable->end()
        if (!iea1->inBondOrAngle(iea2->atom())) {
          bool in14 = iea1->relatedBy14(iea2->atom());
          // 14s are added separately to the Terms
          if (in14 ) {
            if ( !ignore14s ) {
              LOG_ENERGY(BF("Nonbonded interaction between %s - %s in14[%d]\n") % _rep_(iea1->atom()) % _rep_(iea2->atom()) % in14 );
              EnergyNonbond energyNonbond;
              if ( energyNonbond.defineFrom(forceField, in14,
                                            &(*iea1),&(*iea2),this->sharedThis<EnergyNonbond_O>()) )  {
                this->addTerm(energyNonbond);
              }
            }
          } else {
            LOG_ENERGY(BF("Nonbonded interaction between %s - %s in14[%d]\n") % _rep_(iea1->atom()) % _rep_(iea2->atom()) % in14 );
            EnergyNonbond energyNonbond;
            if ( energyNonbond.defineFrom(forceField, in14,
                                          &(*iea1),&(*iea2),this->sharedThis<EnergyNonbond_O>()) )  {
              this->addTerm(energyNonbond);
            }
          }
        }
      }
    }
  } else {
    LOG_ENERGY(BF("There are no non-bonds\n"));
  }
}


void EnergyNonbond_O::constructExcludedAtomListFromAtomTable(AtomTable_sp atomTable, ForceField_sp forceField, bool show_progress)
{
  // ------------------------------------------------------------
  //
  // The new code creates an excluded atom list
  //
  this->_UsesExcludedAtoms = true;
  this->_AtomTable = atomTable;
  this->_FFNonbondDb = forceField->getNonbondDb();
  core::NativeVector_int_mv mv_number_of_excluded_atoms = this->_AtomTable->calculate_excluded_atom_list();
  core::NativeVector_int_sp excluded_atoms_list = mv_number_of_excluded_atoms.second();
  this->_NumberOfExcludedAtomIndices = mv_number_of_excluded_atoms;
  this->_ExcludedAtomIndices = excluded_atoms_list;
}




};
