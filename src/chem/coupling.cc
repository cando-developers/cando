/*
    File: coupling.cc
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


//
// (C) 2004 Christian E. Schafmeister
//



#include <clasp/core/common.h>
#include <cando/chem/coupling.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/monomer.h>
#include <cando/chem/oligomer.h>
#include <clasp/core/package.h>
#include <clasp/core/array.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/plug.h>
#include <cando/chem/bond.h>
//#include	"alchemists.h"
#include <cando/chem/oligomer.h>
#include <clasp/core/wrappers.h>

namespace chem {


void Coupling_O::fields(core::Record_sp node)
{
}

bool	Coupling_O::containsMonomer(Monomer_sp mon)
{
    LOG(BF("checking if coupling(%s) contains monomer: %s") % this->description().c_str() % mon->description().c_str() );
    LOG(BF("About to check monomer1: %s") % this->getMonomer1()->description().c_str() );
    if ( mon == this->getMonomer1() ) { return true; };
    LOG(BF("About to check out monomer: %s") % this->getMonomer2()->description().c_str() );
    if ( mon == this->getMonomer2() ) { return true; };
    LOG(BF("Its neither monomer") );
    return false;
}


#ifdef XML_ARCHIVE
    void	Coupling_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->archiveWeakPointer("oligomer",this->_Oligomer );
    node->attribute("hasError", this->_HasError );
//    node->attribute("statusMessage", this->_StatusMessage);
}
#endif

void DirectionalCoupling_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,source_plug_name), this->_SourcePlugName);
  node->field( INTERN_(kw,target_plug_name), this->_TargetPlugName);
  node->field( INTERN_(kw,source), this->_SourceMonomer);
  node->field( INTERN_(kw,target), this->_TargetMonomer);
}


CL_DEFMETHOD core::T_sp DirectionalCoupling_O::getSourcePlugName() const
{
  return this->_SourcePlugName;
}

CL_DEFMETHOD core::T_sp DirectionalCoupling_O::getTargetPlugName() const
{
  return this->_TargetPlugName;
}

    char DirectionalCoupling_O::otherDirectionalCouplingSide(char s)
{
    if (s == IN_PLUG_PREFIX ) {
        return OUT_PLUG_PREFIX;
    }
    return IN_PLUG_PREFIX;
}

DOCGROUP(cando)
CL_DEFUN bool DirectionalCoupling_O::isCouplingName(core::Symbol_sp name)
{
  core::SimpleString_sp ss = name->symbolName();
  if (ss->length()<1 || !isalpha(ss->rowMajorAref(0).unsafe_character())) return false;
  return true;
}

/*! Strip off the in/out prefix character of a plug name to get a coupling name
    If name is an out-plug-name with digits - the digits become part of the coupling name
    If the name is an in-plug-name and contains a hash - then signal an error 
       - because the coupling name cannot be determined - you will need the coupling to get the coupling name.
 */
DOCGROUP(cando)
CL_DEFUN core::Symbol_sp  DirectionalCoupling_O::couplingName(core::Symbol_sp name)
{
  string sname = name->symbolNameAsString();
  if (sname[0] == IN_PLUG_PREFIX) {
    if (sname.find(IN_PLUG_WILDCARD) != std::string::npos) {
      SIMPLE_ERROR(BF("The coupling name cannot be determined from %s") % _rep_(name));
    }
  } else if (sname[0] != OUT_PLUG_PREFIX) {
    SIMPLE_ERROR(BF("The argument %s must be an in-plug-name or out-plug-name") % _rep_(name));
  }
  std::string rest = sname.substr(1,sname.size());
  return chemkw_intern(rest);
}


DOCGROUP(cando)
CL_DEFUN core::Symbol_sp DirectionalCoupling_O::inPlugName(core::Symbol_sp name)
{
  string sname = name->symbolNameAsString();
  if (!isalpha(sname[0])) {
    SIMPLE_ERROR(BF("You cannot pass a plug name %s to in-plug-name") % _rep_(name));
  }
  size_t dpos = std::string::npos;
  for ( size_t i=0; i<sname.size(); ++i ) {
    if (isdigit(sname[i])) {
      dpos = i;
      break;
    }
  }
  stringstream	ss;
  ss << IN_PLUG_PREFIX;
  if (dpos==std::string::npos) {
    ss << sname;
  } else {
    ss << sname.substr(0,dpos)<< IN_PLUG_WILDCARD;
  }
  return chemkw_intern(ss.str());
};

DOCGROUP(cando)
CL_DEFUN bool DirectionalCoupling_O::outPlugNameMatchesInPlugName(core::Symbol_sp outPlugName, core::Symbol_sp inPlugName) {
  core::SimpleString_sp out = outPlugName->symbolName();
  core::SimpleString_sp in = inPlugName->symbolName();
  if (out->length()<1 || out->rowMajorAref(0).unsafe_character() != OUT_PLUG_PREFIX) {
    SIMPLE_ERROR(BF("The name %s must be an out-plug-name") % _rep_(outPlugName));
  }
  if (in->length()<1 || in->rowMajorAref(0).unsafe_character() != IN_PLUG_PREFIX) {
    SIMPLE_ERROR(BF("The name %s must be an out-plug-name") % _rep_(outPlugName));
  }
  if (out->length() != in->length()) return false;
  size_t pos = 1;
  while (pos<out->length() && pos <in->length()) {
    if (in->rowMajorAref(pos).unsafe_character()==IN_PLUG_WILDCARD) return true;
    if (out->rowMajorAref(pos)!=in->rowMajorAref(pos)) return false;
    pos++;
  }
  return true;
}

DOCGROUP(cando)
CL_DEFUN core::Symbol_sp DirectionalCoupling_O::outPlugName(core::Symbol_sp name)
{
  string sname = name->symbolNameAsString();
  if (!isalpha(sname[0])) {
    SIMPLE_ERROR(BF("You cannot pass a plug name %s to out-plug-name") % _rep_(name));
  }
  if (sname.find(IN_PLUG_WILDCARD)!=std::string::npos) {
    SIMPLE_ERROR(BF("You cannot have a wildcard %s in the input to out-plug-name") % _rep_(name));
  }
  stringstream	ss;
  ss << OUT_PLUG_PREFIX;
  ss << sname;
  return chemkw_intern(ss.str());
};


bool DirectionalCoupling_O::isInPlugName(core::Symbol_sp plugName)
{
  return plugName->symbolName()->rowMajorAref(0).unsafe_character() == IN_PLUG_PREFIX;
};

DOCGROUP(cando)
CL_DEFUN bool chem__is_in_plug_name(core::Symbol_sp plugName)
{
  return DirectionalCoupling_O::isInPlugName(plugName);
}


DOCGROUP(cando)
CL_DEFUN core::Symbol_sp DirectionalCoupling_O::otherPlugName(core::Symbol_sp name)
{
  stringstream	ss;
  ss.str("");
  if ( DirectionalCoupling_O::isInPlugName(name) ) {
    ss << OUT_PLUG_PREFIX;
  } else {
    ss << IN_PLUG_PREFIX;
  }
  ss << DirectionalCoupling_O::couplingName(name)->symbolNameAsString();
  return chemkw_intern(ss.str());
};

//
// Destructor
//

void	DirectionalCoupling_O::checkForErrors(core::T_sp statusMessageStream)
{
  IMPLEMENT_ME();
#if 0
    this->setHasError(false);
    ANN(this->_SourceMonomer);
    if ( this->_SourceMonomer.nilp() )
    {
        this->setHasError(true);
        core::clasp_write_format(BF("SourceMonomer is NULL\n"), statusMessageStream);
    }
    ANN(this->_TargetMonomer);
    if ( this->_TargetMonomer.nilp() )
    {
        this->setHasError(true);
        core::clasp_write_format(BF("TargetMonomer is NULL\n"), statusMessageStream);
    }
#endif
}

CL_LISPIFY_NAME("throwIfBadConnections");
CL_DEFMETHOD void	DirectionalCoupling_O::throwIfBadConnections()
{
  Monomer_sp	mon;
  DirectionalCoupling_sp	coup;
  ANN(this->_SourceMonomer);
  if ( this->_SourceMonomer.nilp() )
  {
    LOG(BF("SourceMonomer is undefined") );
    goto BAD;
  }
  mon = this->_SourceMonomer;
  LOG(BF("SourceMonomer is: %s") % mon->description().c_str()  );
  if ( !mon->hasCouplingWithPlugName(DirectionalCoupling_O::outPlugName(this->getName())))
  {
    LOG(BF("SourceMonomer does not have outPlug with correct name") );
    goto BAD;
  }
  if ( !mon->hasMatchingPlugNameAndCoupling(DirectionalCoupling_O::outPlugName(this->getName()),
                                            this->sharedThis<DirectionalCoupling_O>() ) )
  {
    LOG(BF("SourceMonomers does not have an outPlug that points to us") );
    goto BAD;
  }
  ASSERTNOTNULL(this->_TargetMonomer);
  if ( this->_TargetMonomer.nilp() )
  {
    LOG(BF("TargetMonomer is undefined") );
    goto BAD;
  }
  mon = this->_TargetMonomer;
  LOG(BF("TargetMonomer is: %s") % mon->description().c_str()  );
  if ( !mon->hasCouplingWithPlugName(DirectionalCoupling_O::inPlugName(this->getName())))
  {
    LOG(BF("TargetMonomer does not have inPlug with correct name") );
    goto BAD;
  }
  if ( !mon->hasMatchingPlugNameAndCoupling(DirectionalCoupling_O::inPlugName(this->getName()),
                                            this->sharedThis<DirectionalCoupling_O>() ) )
  {
    LOG(BF("TargetMonomers does not have an inPlug that points to us") );
    goto BAD;
  }
  return;
 BAD:
  SIMPLE_ERROR(BF("Bad coupling %s") % this->sharedThis<DirectionalCoupling_O>()->description() );
}


bool	DirectionalCoupling_O::isInCouplingToMonomer(Monomer_sp mon)
{
		// Ring closing couplings are never in couplings
    ASSERTNOTNULL(this->_TargetMonomer);
    	// If the out Monomer is mon then we are an in coupling to mon
    return this->_TargetMonomer == mon;
};



CL_LISPIFY_NAME("getSourceMonomer");
CL_DEFMETHOD Monomer_sp	DirectionalCoupling_O::getSourceMonomer()
{_OF();
    return this->_SourceMonomer;
};

Monomer_sp	DirectionalCoupling_O::getSourceMonomer_const() const
{_OF();
    return this->_SourceMonomer;
};


CL_LISPIFY_NAME("getTargetMonomer");
CL_DEFMETHOD Monomer_sp	DirectionalCoupling_O::getTargetMonomer()
{_OF();
    return this->_TargetMonomer;
};

Monomer_sp	DirectionalCoupling_O::getTargetMonomer_const() const
{_OF();
    return this->_TargetMonomer;
};



/*!
 * inPlugNames for couplings are the out-plugs for monomers
 */
CL_LISPIFY_NAME("getSourceMonomerPlugName");
CL_DEFMETHOD     core::Symbol_sp	DirectionalCoupling_O::getSourceMonomerPlugName()
{
  return this->_SourcePlugName;
}

CL_LISPIFY_NAME("getName");
CL_DEFMETHOD core::Symbol_sp DirectionalCoupling_O::getName() const {
  return DirectionalCoupling_O::couplingName(this->_SourcePlugName);
}

/*!
 * inPlugNames for couplings are the out-plugs for monomers
 */
CL_LISPIFY_NAME("getTargetMonomerPlugName");
CL_DEFMETHOD     core::Symbol_sp	DirectionalCoupling_O::getTargetMonomerPlugName()
{
  return this->_TargetPlugName;
}



string	DirectionalCoupling_O::description() const
{_OF();
stringstream	ss;
 ss << "DirectionalCoupling( "<< _rep_(this->getName());
    ss << "@";
    ss << std::hex << this << std::dec;
    ss << " ";
    ANN(this->_SourceMonomer);
    if ( this->_SourceMonomer.nilp())
    {
	ss << "-NULL_Monomer-";
    } else
    {
	ss << this->getSourceMonomer_const()->description();
    }
    ss << "->";
    ANN(this->_TargetMonomer);
    if ( this->_TargetMonomer.nilp() )
    {
	ss << "-NULL_Monomer-";
    } else
    {
	ss << this->getTargetMonomer_const()->description();
    }
    ss << ")";
    return ss.str();
}



CL_LISPIFY_NAME("setPlugNames");
CL_DEFMETHOD     void	DirectionalCoupling_O::setPlugNames(core::Symbol_sp sourcePlugName, core::Symbol_sp targetPlugName)
{
  this->_SourcePlugName = sourcePlugName;
  this->_TargetPlugName = targetPlugName;
  if ( this->_SourceMonomer.notnilp() )
  {
    this->_SourceMonomer->fixPlugNameForCoupling(this->sharedThis<DirectionalCoupling_O>());
  }
  if ( this->_TargetMonomer.notnilp() )
  {
    this->_TargetMonomer->fixPlugNameForCoupling(this->sharedThis<DirectionalCoupling_O>());
  }
}


void	DirectionalCoupling_O::setSourceMonomer_NoSignal(Monomer_sp m)
{
    ASSERTNOTNULLP(m,"above monomer is undefined");
    this->_SourceMonomer = m->sharedThis<Monomer_O>();
    LOG(BF("Setting in monomer to (%s)") % m->description().c_str()  );
}

void	DirectionalCoupling_O::setTargetMonomer_NoSignal(Monomer_sp m)
{
    ASSERTNOTNULLP(m,"below monomer is undefined");
    this->_TargetMonomer = m->sharedThis<Monomer_O>();
    LOG(BF("Setting out monomer to (%s)") % m->description().c_str()  );
}



CL_LISPIFY_NAME("setSourceMonomer");
CL_DEFMETHOD void	DirectionalCoupling_O::setSourceMonomer(Monomer_sp m)
{
    this->setSourceMonomer_NoSignal(m);
//    SIGNAL(this,Coupling_connectionsChanged);
}

CL_LISPIFY_NAME("setTargetMonomer");
CL_DEFMETHOD void DirectionalCoupling_O::setTargetMonomer(Monomer_sp m)
{
    this->setTargetMonomer_NoSignal(m);
//    SIGNAL(this,Coupling_connectionsChanged);
}




DirectionalCouplingSide DirectionalCoupling_O::couplingSideOfMonomer( Monomer_sp mon )
{_OF();
    ASSERTNOTNULL( mon );
    ASSERT(mon.notnilp());
    if ( mon == this->getSourceMonomer() ) return InSide;
    if ( mon == this->getTargetMonomer() ) return OutSide;
    return NoSide;
}





CL_LISPIFY_NAME("getOtherSideMonomer");
CL_DEFMETHOD Monomer_sp	DirectionalCoupling_O::getOtherSideMonomer(Monomer_sp mon)
{
    LOG(BF("About to check in monomer") );
    if ( mon == this->getSourceMonomer() ) { return this->getTargetMonomer(); };
    LOG(BF("About to check out monomer") );
    if ( mon == this->getTargetMonomer() ) { return this->getSourceMonomer(); };
    LOG(BF("Its neither monomer") );
    SIMPLE_ERROR(BF("Monomer is not in coupling %s") % mon->description());
}






void	DirectionalCoupling_O::doCoupling(Residue_sp inResidue, Residue_sp outResidue )
{
Monomer_sp	min, mout;
PlugWithMates_sp	minPlug, moutPlug;
Residue_sp	weakInRes, weakOutRes;
Bond_sp		b;
Atom_sp		inB0, outB0, inB1, outB1;
    min = this->getSourceMonomer();
    mout = this->getTargetMonomer();
    minPlug = min->getPlugNamed(this->getSourceMonomerPlugName()).as<PlugWithMates_O>();
    moutPlug = mout->getPlugNamed(this->getTargetMonomerPlugName()).as<PlugWithMates_O>();
    ASSERTNOTNULL(minPlug);
    ASSERTNOTNULL(moutPlug);
    LOG(BF("Everything seems ok with the residues, let's do the coupling") );
    	// connect the atoms in min/mout according to minPlug/moutPlug
    inB0 = gc::As_unsafe<Atom_sp>(inResidue->atomWithName(minPlug->getB0()));
    outB0 = gc::As_unsafe<Atom_sp>(outResidue->atomWithName(moutPlug->getB0()));
    b = inB0->bondTo(outB0,singleBond);
    if ( moutPlug->getB1().notnilp() ) {
      outB1 = gc::As_unsafe<Atom_sp>(outResidue->atomWithName(moutPlug->getB1()));
	if ( minPlug->getB1().notnilp() ) {
          inB1 = gc::As_unsafe<Atom_sp>(inResidue->atomWithName(minPlug->getB1()));
	 	//
		// Create the second bond
		//
	    b = inB1->bondTo(outB1,singleBond);
	} else {
		//
		// We have in(B0) and out(B0)/out(B1)
		// Connect both out atoms back to B0
		//
	    b = inB0->bondTo(outB1,singleBond);
	}
    } else if ( minPlug->getB1().notnilp() ) {
		//
		// Now we have in(B0)/in(B1) and only out(B0)
		//
	b = inB1->bondTo(outB0,singleBond);
    }
    return;
}


Plug_sp	DirectionalCoupling_O::getPlugForMonomer(Monomer_sp mon)
{
    if ( this->_SourceMonomer == mon ) {
	core::Symbol_sp plugName = this->getSourceMonomerPlugName();
	ASSERT(this->_SourceMonomer.notnilp());
	return this->_SourceMonomer->getPlugNamed(plugName);
    }
    core::Symbol_sp plugName = this->getTargetMonomerPlugName();
    ASSERT(this->_TargetMonomer.notnilp());
    return this->_TargetMonomer->getPlugNamed(plugName);
}

Plug_sp	DirectionalCoupling_O::getPlugForOtherMonomer(Monomer_sp mon)
{
    if ( this->_SourceMonomer == mon ) {
	core::Symbol_sp plugName = this->getTargetMonomerPlugName();
	ASSERT(this->_TargetMonomer.notnilp());
	return this->_TargetMonomer->getPlugNamed(plugName);
    }
    core::Symbol_sp plugName = this->getSourceMonomerPlugName();
    ASSERT(this->_SourceMonomer.notnilp());
    return this->_SourceMonomer->getPlugNamed(plugName);
}







void	DirectionalCoupling_O::couple( Monomer_sp sin, Monomer_sp sout )
{_OF();
DirectionalCoupling_sp	me;
SIMPLE_ERROR(BF("I wasn't sure if this was ever called.  Take out this THROW if it is"));
    this->setSourceMonomer_NoSignal(sin);
    this->setTargetMonomer_NoSignal(sout);
    me = this->sharedThis<DirectionalCoupling_O>();
    core::Symbol_sp inMonomerPlugName = DirectionalCoupling_O::outPlugName(this->getName());
    core::Symbol_sp outMonomerPlugName = DirectionalCoupling_O::inPlugName(this->getName());
    sin->addCoupling(inMonomerPlugName,me);
    sout->addCoupling(outMonomerPlugName,me);
//    SIGNAL(this,Coupling_connectionsChanged);
};

CL_NAME(chem:|Coupling_O::isRingClosing|);
CL_DEFMETHOD bool Coupling_O::isRingClosing() const
{
  return false;
}


Coupling_sp DirectionalCoupling_O::copyAndReplaceMonomers( Monomer_sp monomer1, Monomer_sp monomer2 ) {
  return DirectionalCoupling_O::make( monomer1, this->_SourcePlugName, this->_TargetPlugName, monomer2 );
}

#ifdef XML_ARCHIVE
void RingCoupling_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("plug1", this->_Plug1);
    node->attribute("plug2", this->_Plug2);
    node->archiveWeakPointer("monomer1", this->_Monomer1 );
    node->archiveWeakPointer("monomer2", this->_Monomer2 );
}
#endif

void RingCoupling_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,plug1), this->_Plug1);
  node->field( INTERN_(kw,monomer1), this->_Monomer1);
  node->field( INTERN_(kw,plug2), this->_Plug2);
  node->field( INTERN_(kw,monomer2), this->_Monomer2);
}


    core::Symbol_sp RingCoupling_O::getName() const
{
    string name = "ringCoupling["+this->_Plug1->symbolNameAsString()+">"+this->_Plug2->symbolNameAsString()+"]";
    return chemkw_intern(name);
}


CL_LISPIFY_NAME("setMonomer1");
CL_DEFMETHOD void	RingCoupling_O::setMonomer1(Monomer_sp mon)
{
    this->_Monomer1 = mon;
}

CL_LISPIFY_NAME("getMonomer1");
CL_DEFMETHOD Monomer_sp RingCoupling_O::getMonomer1()
{_OF();
    ASSERTNOTNULLP(this->_Monomer1,
    	"SourceMonomer is NULL for coupling: "+this->description());
    return this->_Monomer1;
}

Monomer_sp RingCoupling_O::getMonomer1_const() const
{_OF();
    ASSERTNOTNULLP(this->_Monomer1,
    	"SourceMonomer is NULL for coupling: "+this->description());
    return this->_Monomer1;
}


CL_LISPIFY_NAME("setMonomer2");
CL_DEFMETHOD void	RingCoupling_O::setMonomer2(Monomer_sp mon)
{
    this->_Monomer2 = mon;
}

CL_LISPIFY_NAME("getMonomer2");
CL_DEFMETHOD Monomer_sp RingCoupling_O::getMonomer2()
{_OF();
    ASSERTNOTNULLP(this->_Monomer2,
    	"Monomer2 is NULL for coupling: "+this->description());
    return this->_Monomer2;
}

Monomer_sp RingCoupling_O::getMonomer2_const() const
{_OF();
    ASSERTNOTNULLP(this->_Monomer2,
    	"Monomer2 is NULL for coupling: "+this->description());
    return this->_Monomer2;
}

CL_LISPIFY_NAME("getOtherSideMonomer");
CL_DEFMETHOD Monomer_sp	RingCoupling_O::getOtherSideMonomer(Monomer_sp mon)
{
    LOG(BF("About to check in monomer") );
    if ( mon == this->getMonomer1() ) { return this->getMonomer2(); };
    LOG(BF("About to check out monomer") );
    if ( mon == this->getMonomer2() ) { return this->getMonomer1(); };
    LOG(BF("Its neither monomer") );
    SIMPLE_ERROR(BF("Monomer is not in coupling %s" ) % mon->description());
}



void	RingCoupling_O::checkForErrors(core::T_sp errorStream)
{
  IMPLEMENT_ME();
#if 0
  ANN(this->_Monomer1);
    if ( this->_Monomer1.nilp() )
    {
        this->setHasError(true);
        core::clasp_write_format(BF("Monomer1 is NULL\n"), errorStream);
    }
    ANN(this->_Monomer2);
    if ( this->_Monomer2.nilp() )
    {
        this->setHasError(true);
        core::clasp_write_format(BF("Monomer2 is NULL\n"), errorStream);
    }
#endif
}

void	RingCoupling_O::throwIfBadConnections()
{
Monomer_sp	mon;
RingCoupling_sp	coup;
    ANN(this->_Monomer1);
    if ( this->_Monomer1.nilp() )
    {
        LOG(BF("Monomer1 is undefined") );
        goto BAD;
    }
    mon = this->_Monomer1;
    LOG(BF("Monomer1 is: %s") % mon->description().c_str()  );
    if ( !mon->hasCouplingWithPlugName(DirectionalCoupling_O::outPlugName(this->getPlug1())))
    {
        LOG(BF("Monomer1 does not have outPlug with correct name") );
        goto BAD;
    }
    if ( !mon->hasMatchingPlugNameAndCoupling(DirectionalCoupling_O::outPlugName(this->getPlug1()),
    					this->sharedThis<RingCoupling_O>() ) )
    {
        LOG(BF("Monomer1s does not have an outPlug that points to us") );
        goto BAD;
    }
    ASSERTNOTNULL(this->_Monomer2);
    if ( this->_Monomer2.nilp() )
    {
        LOG(BF("Monomer2 is undefined") );
        goto BAD;
    }
    mon = this->_Monomer2;
    LOG(BF("Monomer2 is: %s") % mon->description().c_str()  );
    if ( !mon->hasCouplingWithPlugName(DirectionalCoupling_O::outPlugName(this->getPlug2())))
    {
        LOG(BF("Monomer2 does not have outPlug with correct name") );
        goto BAD;
    }
    if ( !mon->hasMatchingPlugNameAndCoupling(DirectionalCoupling_O::outPlugName(this->getPlug2()),
    					this->sharedThis<RingCoupling_O>() ) )
    {
        LOG(BF("Monomer2s does not have an inPlug that points to us") );
        goto BAD;
    }
    return;
BAD:
    SIMPLE_ERROR(BF("Bad coupling %s" ) % this->sharedThis<RingCoupling_O>()->description() );
}


void RingCoupling_O::setPlug1(core::Symbol_sp p)
{
    this->_Plug1 = p;
}

void RingCoupling_O::setPlug2(core::Symbol_sp p)
{
    this->_Plug2 = p;
}

CL_LISPIFY_NAME("getPlug1");
CL_DEFMETHOD core::Symbol_sp RingCoupling_O::getPlug1()
{_OF();
    return this->_Plug1;
}
CL_LISPIFY_NAME("getPlug2");
CL_DEFMETHOD core::Symbol_sp RingCoupling_O::getPlug2()
{_OF();
    return this->_Plug2;
}

string	RingCoupling_O::description() const
{_OF();
stringstream	ss;
 ss << "RingCoupling( "<< _rep_(this->getName());
    ss << "@";
    ss << std::hex << this << std::dec;
    ss << " ";
    ANN(this->_Monomer1);
    if ( this->_Monomer1.nilp())
    {
	ss << "-NULL_Monomer-";
    } else
    {
	ss << this->getMonomer1_const()->description();
    }
    ss << "->";
    ANN(this->_Monomer2);
    if ( this->_Monomer2.nilp() )
    {
	ss << "-NULL_Monomer-";
    } else
    {
	ss << this->getMonomer2_const()->description();
    }
    ss << ")";
    return ss.str();
}



void	RingCoupling_O::doCoupling(Residue_sp inResidue, Residue_sp outResidue )
{
Monomer_sp	mon1, mon2;
PlugWithMates_sp	mon1Plug, mon2Plug;
Residue_sp	weakInRes, weakOutRes;
Bond_sp		b;
Atom_sp		inB0, outB0, inB1, outB1;
    mon1 = this->getMonomer1();
    mon2 = this->getMonomer2();
    mon1Plug = mon1->getPlugNamed(this->getPlug1()).as<PlugWithMates_O>();
    mon2Plug = mon2->getPlugNamed(this->getPlug2()).as<PlugWithMates_O>();
    ASSERTNOTNULL(mon1Plug);
    ASSERTNOTNULL(mon2Plug);
    LOG(BF("Everything seems ok with the residues, let's do the coupling") );
    	// connect the atoms in mon1/mon2 according to mon1Plug/mon2Plug
    inB0 = gc::As_unsafe<Atom_sp>(inResidue->atomWithName(mon1Plug->getB0()));
    outB0 = gc::As_unsafe<Atom_sp>(outResidue->atomWithName(mon2Plug->getB0()));
    b = inB0->bondTo(outB0,singleBond);
    if ( mon2Plug->getB1().notnilp() ) {
      outB1 = gc::As_unsafe<Atom_sp>(outResidue->atomWithName(mon2Plug->getB1()));
	if ( mon1Plug->getB1().notnilp() ) {
          inB1 = gc::As_unsafe<Atom_sp>(inResidue->atomWithName(mon1Plug->getB1()));
	 	//
		// Create the second bond
		//
	    b = inB1->bondTo(outB1,singleBond);
	} else {
		//
		// We have in(B0) and out(B0)/out(B1)
		// Connect both out atoms back to B0
		//
	    b = inB0->bondTo(outB1,singleBond);
	}
    } else if ( mon1Plug->getB1().notnilp() ) {
		//
		// Now we have in(B0)/in(B1) and only out(B0)
		//
	b = inB1->bondTo(outB0,singleBond);
    }
    return;
}


Coupling_sp RingCoupling_O::copyAndReplaceMonomers( Monomer_sp monomer1, Monomer_sp monomer2 ) {
  return RingCoupling_O::make( monomer1, this->_Plug1, monomer2, this->_Plug2 );
}



















};
