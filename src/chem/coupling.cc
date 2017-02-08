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
#include <cando/chem/candoDatabase.h>
#include <cando/chem/plug.h>
#include <cando/chem/bond.h>
//#include	"alchemists.h"
#include <cando/chem/oligomer.h>
#include <clasp/core/wrappers.h>

namespace chem {




void Coupling_O::initialize()
{
    this->Base::initialize();
    this->_Oligomer = _Nil<core::T_O>();
    this->_HasError = false;
    this->_Selected = false;
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


void	Coupling_O::setOligomer(Oligomer_sp o)
{
    this->_Oligomer = o->sharedThis<Oligomer_O>();
};


void	Coupling_O::resetOligomer()
{
  this->_Oligomer = _Nil<core::T_O>();
}

Oligomer_sp	Coupling_O::getOligomer()
{_OF();
    ASSERTNOTNULL(this->_Oligomer);
    return this->_Oligomer;
};








    char DirectionalCoupling_O::otherDirectionalCouplingSide(char s)
{
    if (s == IN_PLUG_PREFIX ) {
        return OUT_PLUG_PREFIX;
    }
    return IN_PLUG_PREFIX;
}

/*! Strip off the in/out prefix character of a plug name to get a coupling name
 */
CL_DEFUN core::Symbol_sp  DirectionalCoupling_O::couplingName(core::Symbol_sp name)
    {
	string sname = name->symbolNameAsString();
	LOG(BF("couplingName symbol[%s]  symbolNameAsString[%s]") % _rep_(name) % sname );
	string	rest;
	if ( sname[0] == IN_PLUG_PREFIX || sname[0] == OUT_PLUG_PREFIX ) {
	    rest = sname.substr(1,sname.size());
	} else rest = sname;
	return chemkw_intern(rest);
    }

CL_DEFUN core::Symbol_sp DirectionalCoupling_O::inPlugName(core::Symbol_sp name)
    {
	stringstream	ss;
	ss.str("");
	ss << IN_PLUG_PREFIX;
	ss << DirectionalCoupling_O::couplingName(name)->symbolNameAsString();
	return chemkw_intern(ss.str());
    };


CL_DEFUN core::Symbol_sp DirectionalCoupling_O::outPlugName(core::Symbol_sp name)
{
stringstream	ss;
    ss.str("");
    ss << OUT_PLUG_PREFIX;
    ss << DirectionalCoupling_O::couplingName(name)->symbolNameAsString();
    return chemkw_intern(ss.str());
};


bool DirectionalCoupling_O::isInPlugName(core::Symbol_sp plugName)
{
  return plugName->symbolName()->rowMajorAref(0).unsafe_character() == IN_PLUG_PREFIX;
};

CL_DEFUN bool chem__is_in_plug_name(core::Symbol_sp plugName)
{
  return DirectionalCoupling_O::isInPlugName(plugName);
}


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











void DirectionalCoupling_O::initialize()
{
    this->Base::initialize();
    this->_InMonomer = _Nil<core::T_O>();
    this->_OutMonomer = _Nil<core::T_O>();
    this->_Name = _Nil<core::Symbol_O>();
}

//
// Destructor
//


void	DirectionalCoupling_O::resetIn()
{
  this->_InMonomer = _Nil<core::T_O>();
};


void	DirectionalCoupling_O::resetOut()
{
  this->_OutMonomer = _Nil<core::T_O>();
};




    void	DirectionalCoupling_O::checkForErrors(core::T_sp statusMessageStream)
{
    this->setHasError(false);
    ANN(this->_InMonomer);
    if ( this->_InMonomer.nilp() )
    {
        this->setHasError(true);
        core::clasp_write_format(BF("InMonomer is NULL\n"), statusMessageStream);
    }
    ANN(this->_OutMonomer);
    if ( this->_OutMonomer.nilp() )
    {
        this->setHasError(true);
        core::clasp_write_format(BF("OutMonomer is NULL\n"), statusMessageStream);
    }
}

CL_LISPIFY_NAME("throwIfBadConnections");
CL_DEFMETHOD void	DirectionalCoupling_O::throwIfBadConnections()
{
Monomer_sp	mon;
DirectionalCoupling_sp	coup;
    ANN(this->_InMonomer);
    if ( this->_InMonomer.nilp() )
    {
        LOG(BF("InMonomer is undefined") );
        goto BAD;
    }
    mon = this->_InMonomer;
    LOG(BF("InMonomer is: %s") % mon->description().c_str()  );
    if ( !mon->hasCouplingWithPlugName(DirectionalCoupling_O::outPlugName(this->getName())))
    {
        LOG(BF("InMonomer does not have outPlug with correct name") );
        goto BAD;
    }
    if ( !mon->hasMatchingPlugNameAndCoupling(DirectionalCoupling_O::outPlugName(this->getName()),
    					this->sharedThis<DirectionalCoupling_O>() ) )
    {
        LOG(BF("InMonomers does not have an outPlug that points to us") );
        goto BAD;
    }
    ASSERTNOTNULL(this->_OutMonomer);
    if ( this->_OutMonomer.nilp() )
    {
        LOG(BF("OutMonomer is undefined") );
        goto BAD;
    }
    mon = this->_OutMonomer;
    LOG(BF("OutMonomer is: %s") % mon->description().c_str()  );
    if ( !mon->hasCouplingWithPlugName(DirectionalCoupling_O::inPlugName(this->getName())))
    {
        LOG(BF("OutMonomer does not have inPlug with correct name") );
        goto BAD;
    }
    if ( !mon->hasMatchingPlugNameAndCoupling(DirectionalCoupling_O::inPlugName(this->getName()),
    					this->sharedThis<DirectionalCoupling_O>() ) )
    {
        LOG(BF("OutMonomers does not have an inPlug that points to us") );
        goto BAD;
    }
    return;
BAD:
    SIMPLE_ERROR(BF("Bad coupling %s") % this->sharedThis<DirectionalCoupling_O>()->description() );
}


#ifdef XML_ARCHIVE
    void	DirectionalCoupling_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("name", this->_Name );
    node->archiveWeakPointer("inMonomer", this->_InMonomer );
    node->archiveWeakPointer("outMonomer", this->_OutMonomer );
}
#endif


bool	DirectionalCoupling_O::isInCouplingToMonomer(Monomer_sp mon)
{
		// Ring closing couplings are never in couplings
    ASSERTNOTNULL(this->_OutMonomer);
    	// If the out Monomer is mon then we are an in coupling to mon
    return this->_OutMonomer == mon;
};



CL_LISPIFY_NAME("getInMonomer");
CL_DEFMETHOD Monomer_sp	DirectionalCoupling_O::getInMonomer()
{_OF();
    ASSERTNOTNULLP(this->_InMonomer,
    	"InMonomer is NULL for coupling: "+this->description());
    return this->_InMonomer;
};

Monomer_sp	DirectionalCoupling_O::getInMonomer_const() const
{_OF();
    ASSERTNOTNULLP(this->_InMonomer,
    	"InMonomer is NULL for coupling: "+this->description());
    return this->_InMonomer;
};


CL_LISPIFY_NAME("getOutMonomer");
CL_DEFMETHOD Monomer_sp	DirectionalCoupling_O::getOutMonomer()
{_OF();
    ASSERTNOTNULLP(this->_OutMonomer,
    	"OutMonomer is NULL for coupling: "+this->description());
    return this->_OutMonomer;
};

Monomer_sp	DirectionalCoupling_O::getOutMonomer_const() const
{_OF();
    ASSERTNOTNULLP(this->_OutMonomer,
    	"OutMonomer is NULL for coupling: "+this->description());
    return this->_OutMonomer;
};



/*!
 * inPlugNames for couplings are the out-plugs for monomers
 */
CL_LISPIFY_NAME("getInMonomerPlugName");
CL_DEFMETHOD     core::Symbol_sp	DirectionalCoupling_O::getInMonomerPlugName()
{
    return DirectionalCoupling_O::outPlugName(this->_Name);
}


/*!
 * inPlugNames for couplings are the out-plugs for monomers
 */
CL_LISPIFY_NAME("getOutMonomerPlugName");
CL_DEFMETHOD     core::Symbol_sp	DirectionalCoupling_O::getOutMonomerPlugName()
{
    return DirectionalCoupling_O::inPlugName(this->_Name);
}



string	DirectionalCoupling_O::description() const
{_OF();
stringstream	ss;
 ss << "DirectionalCoupling( "<< _rep_(this->getName());
    ss << "@";
    ss << std::hex << this << std::dec;
    ss << " ";
    ANN(this->_InMonomer);
    if ( this->_InMonomer.nilp())
    {
	ss << "-NULL_Monomer-";
    } else
    {
	ss << this->getInMonomer_const()->description();
    }
    ss << "->";
    ANN(this->_OutMonomer);
    if ( this->_OutMonomer.nilp() )
    {
	ss << "-NULL_Monomer-";
    } else
    {
	ss << this->getOutMonomer_const()->description();
    }
    ss << ")";
    return ss.str();
}




CL_LISPIFY_NAME("setName");
CL_DEFMETHOD     void	DirectionalCoupling_O::setName(core::Symbol_sp nm)
{
    this->_Name = DirectionalCoupling_O::couplingName(nm);
    ANN(this->_InMonomer);
    if ( this->_InMonomer.notnilp() )
    {
	this->_InMonomer->fixPlugNameForCoupling(this->sharedThis<DirectionalCoupling_O>());
    }
    ANN(this->_OutMonomer);
    if ( this->_OutMonomer.notnilp() )
    {
	this->_OutMonomer->fixPlugNameForCoupling(this->sharedThis<DirectionalCoupling_O>());
    }
}


void	DirectionalCoupling_O::setInMonomer_NoSignal(Monomer_sp m)
{
    ASSERTNOTNULLP(m,"above monomer is undefined");
    this->_InMonomer = m->sharedThis<Monomer_O>();
    LOG(BF("Setting in monomer to (%s)") % m->description().c_str()  );
}

void	DirectionalCoupling_O::setOutMonomer_NoSignal(Monomer_sp m)
{
    ASSERTNOTNULLP(m,"below monomer is undefined");
    this->_OutMonomer = m->sharedThis<Monomer_O>();
    LOG(BF("Setting out monomer to (%s)") % m->description().c_str()  );
}



CL_LISPIFY_NAME("setInMonomer");
CL_DEFMETHOD void	DirectionalCoupling_O::setInMonomer(Monomer_sp m)
{
    this->setInMonomer_NoSignal(m);
//    SIGNAL(this,Coupling_connectionsChanged);
}

CL_LISPIFY_NAME("setOutMonomer");
CL_DEFMETHOD void	DirectionalCoupling_O::setOutMonomer(Monomer_sp m)
{
    this->setOutMonomer_NoSignal(m);
//    SIGNAL(this,Coupling_connectionsChanged);
}




DirectionalCouplingSide DirectionalCoupling_O::couplingSideOfMonomer( Monomer_sp mon )
{_OF();
    ASSERTNOTNULL( mon );
    ASSERT(mon.notnilp());
    if ( mon == this->getInMonomer() ) return InSide;
    if ( mon == this->getOutMonomer() ) return OutSide;
    return NoSide;
}





CL_LISPIFY_NAME("getOtherSideMonomer");
CL_DEFMETHOD Monomer_sp	DirectionalCoupling_O::getOtherSideMonomer(Monomer_sp mon)
{
    LOG(BF("About to check in monomer") );
    if ( mon == this->getInMonomer() ) { return this->getOutMonomer(); };
    LOG(BF("About to check out monomer") );
    if ( mon == this->getOutMonomer() ) { return this->getInMonomer(); };
    LOG(BF("Its neither monomer") );
    SIMPLE_ERROR(BF("Monomer is not in coupling %s") % mon->description());
}






void	DirectionalCoupling_O::doCoupling(Residue_sp inResidue, Residue_sp outResidue )
{
Monomer_sp	min, mout;
PlugWithMates_sp	minPlug, moutPlug;
Residue_wp	weakInRes, weakOutRes;
Bond_sp		b;
Atom_sp		inB0, outB0, inB1, outB1;
    if ( this->getHasError() ) return;
    min = this->getInMonomer();
    mout = this->getOutMonomer();
    if ( min->getHasError() ) return;
    if ( mout->getHasError() ) return;
    ASSERT(min.notnilp());
    ASSERT(mout.notnilp());
    minPlug = min->getPlugNamed(this->getInMonomerPlugName()).as<PlugWithMates_O>();
    moutPlug = mout->getPlugNamed(this->getOutMonomerPlugName()).as<PlugWithMates_O>();
    ASSERTNOTNULL(minPlug);
    ASSERTNOTNULL(moutPlug);
    		//
		// Each monomer has a temporary residue must be identical to the
		// residues that are passed to this routine.
		// This is to ensure that the residues that we are doing the coupling
		// between are the residues that represent the monomers.
		//
    weakInRes = min->getTemporaryResidue();
    ASSERTNOTNULL(weakInRes);
    	// make sure the monomers residue is the same as the one the caller passes us
    ASSERT(weakInRes == inResidue );
    weakOutRes = mout->getTemporaryResidue();
    ASSERTNOTNULL(weakOutRes);
    	// make sure the monomers residue is the same as the one the caller passes us
    ASSERT(weakOutRes == outResidue );
    LOG(BF("Everything seems ok with the residues, let's do the coupling") );
    	// connect the atoms in min/mout according to minPlug/moutPlug
    inB0 = inResidue->atomWithName(minPlug->getB0());
    outB0 = outResidue->atomWithName(moutPlug->getB0());
    b = inB0->bondTo(outB0,singleBond);
    if ( moutPlug->getB1().notnilp() ) {
	outB1 = outResidue->atomWithName(moutPlug->getB1());
	if ( minPlug->getB1().notnilp() ) {
	    inB1 = inResidue->atomWithName(minPlug->getB1());
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
    if ( this->_InMonomer == mon ) {
	core::Symbol_sp plugName = this->getInMonomerPlugName();
	ASSERT(this->_InMonomer.notnilp());
	return this->_InMonomer->getPlugNamed(plugName);
    }
    core::Symbol_sp plugName = this->getOutMonomerPlugName();
    ASSERT(this->_OutMonomer.notnilp());
    return this->_OutMonomer->getPlugNamed(plugName);
}

Plug_sp	DirectionalCoupling_O::getPlugForOtherMonomer(Monomer_sp mon)
{
    if ( this->_InMonomer == mon ) {
	core::Symbol_sp plugName = this->getOutMonomerPlugName();
	ASSERT(this->_OutMonomer.notnilp());
	return this->_OutMonomer->getPlugNamed(plugName);
    }
    core::Symbol_sp plugName = this->getInMonomerPlugName();
    ASSERT(this->_InMonomer.notnilp());
    return this->_InMonomer->getPlugNamed(plugName);
}







void	DirectionalCoupling_O::couple( Monomer_sp sin, Monomer_sp sout )
{_OF();
DirectionalCoupling_sp	me;
SIMPLE_ERROR(BF("I wasn't sure if this was ever called.  Take out this THROW if it is"));
    this->setInMonomer_NoSignal(sin);
    this->setOutMonomer_NoSignal(sout);
    me = this->sharedThis<DirectionalCoupling_O>();
    core::Symbol_sp inMonomerPlugName = DirectionalCoupling_O::outPlugName(this->getName());
    core::Symbol_sp outMonomerPlugName = DirectionalCoupling_O::inPlugName(this->getName());
    sin->addCoupling(inMonomerPlugName,me);
    sout->addCoupling(outMonomerPlugName,me);
//    SIGNAL(this,Coupling_connectionsChanged);
};







void RingCoupling_O::initialize()
{_OF();
    this->Base::initialize();
    this->_Monomer1 = _Nil<core::T_O>();
    this->_Monomer2 = _Nil<core::T_O>();
    this->_Plug1 = _Nil<core::Symbol_O>();
    this->_Plug2 = _Nil<core::Symbol_O>();
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
    	"InMonomer is NULL for coupling: "+this->description());
    return this->_Monomer1;
}

Monomer_sp RingCoupling_O::getMonomer1_const() const
{_OF();
    ASSERTNOTNULLP(this->_Monomer1,
    	"InMonomer is NULL for coupling: "+this->description());
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
    this->setHasError(false);
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
    return DirectionalCoupling_O::outPlugName(this->_Plug1);
}
CL_LISPIFY_NAME("getPlug2");
CL_DEFMETHOD core::Symbol_sp RingCoupling_O::getPlug2()
{_OF();
    return DirectionalCoupling_O::outPlugName(this->_Plug2);
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
Residue_wp	weakInRes, weakOutRes;
Bond_sp		b;
Atom_sp		inB0, outB0, inB1, outB1;
    if ( this->getHasError() ) return;
    mon1 = this->getMonomer1();
    mon2 = this->getMonomer2();
    if ( mon1->getHasError() ) return;
    if ( mon2->getHasError() ) return;
    ASSERT(mon1.notnilp());
    ASSERT(mon2.notnilp());
    mon1Plug = mon1->getPlugNamed(this->getPlug1()).as<PlugWithMates_O>();
    mon2Plug = mon2->getPlugNamed(this->getPlug2()).as<PlugWithMates_O>();
    ASSERTNOTNULL(mon1Plug);
    ASSERTNOTNULL(mon2Plug);
    		//
		// Each monomer has a temporary residue must be identical to the
		// residues that are passed to this routine.
		// This is to ensure that the residues that we are doing the coupling
		// between are the residues that represent the monomers.
		//
    weakInRes = mon1->getTemporaryResidue();
    ASSERTNOTNULL(weakInRes);
    	// make sure the monomers residue is the same as the one the caller passes us
    ASSERT(weakInRes == inResidue );
    weakOutRes = mon2->getTemporaryResidue();
    ASSERTNOTNULL(weakOutRes);
    	// make sure the monomers residue is the same as the one the caller passes us
    ASSERT(weakOutRes == outResidue );
    LOG(BF("Everything seems ok with the residues, let's do the coupling") );
    	// connect the atoms in mon1/mon2 according to mon1Plug/mon2Plug
    inB0 = inResidue->atomWithName(mon1Plug->getB0());
    outB0 = outResidue->atomWithName(mon2Plug->getB0());
    b = inB0->bondTo(outB0,singleBond);
    if ( mon2Plug->getB1().notnilp() ) {
	outB1 = outResidue->atomWithName(mon2Plug->getB1());
	if ( mon1Plug->getB1().notnilp() ) {
	    inB1 = inResidue->atomWithName(mon1Plug->getB1());
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





















};
