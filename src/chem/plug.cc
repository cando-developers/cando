/*
    File: plug.cc
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

#include <clasp/core/common.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/plug.h>
#include <cando/chem/candoDatabase.h>
#include <clasp/core/cons.h>
#include <cando/chem/coupling.h>
#include <clasp/core/environment.h>
#include <clasp/core/numbers.h>
#include <cando/chem/constitution.h>
#include <clasp/core/wrappers.h>


namespace chem
{



    core::Symbol_sp RingClosingMate_O::getName() const
    {
	return _lisp->intern("UnnamedRingClosingMate");
    };


//
// Constructor
//




    string RingClosingMate_O::descriptionOfContents() const
    {
	stringstream ss, shex;
	return ss.str();
    }

#ifdef XML_ARCHIVE
    void	RingClosingMate_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
    }
#endif

/*
  __BEGIN_DOC(classes.RingClosingMate)
  \inheritsFrom{EntityNameSetWithCap}
  A EntityNameSet that keeps track of a capping monomer that is used to cap training oligomers when they are being defined.
  The capping monomer is supposed to be small and best represent the other members of the RingClosingMate.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.RingClosingMate.!class.RingClosingMate)
  \requiredKeyed{groupNames:}{List::groupNames}

  Initialize a RingClosingMate object.
  __END_DOC
*/

    core::Symbol_sp Mate_O::getName() const
    {
	return _Nil<core::Symbol_O>();
    }


#if 0
    Mate_sp Mate_O::create(CandoDatabase_sp db)
    {
	GC_ALLOCATE(Mate_O, mate );
	mate->setCandoDatabase(db);
	return mate;
    }
#endif

//
// Constructor
//



    string Mate_O::descriptionOfContents() const
    {
	stringstream ss, shex;
	ss << "( Cap "<< this->_Cap<<") ";
	return ss.str();
    }

#ifdef XML_ARCHIVE
    void	Mate_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute( "cap", this->_Cap );
    }
#endif

/*
  __BEGIN_DOC(classes.Mate)
  \inheritsFrom{EntityNameSetWithCap}
  A EntityNameSet that keeps track of a capping monomer that is used to cap training oligomers when they are being defined.
  The capping monomer is supposed to be small and best represent the other members of the Mate.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.Mate.!class.Mate)
  \requiredKeyed{cap:}{Text::capName}
  \requiredKeyed{groupNames:}{List::groupNames}

  Initialize a Mate object.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.Plug)
  Defines one or two atoms of this monomer that can be plugged into, a plug name and a
  collection of Mate objects that can act as mates for this plug.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.Plug.!class.Plug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a Plug object. Plugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  Outgoing plugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  __END_DOC
*/


    string	Plug_O::descriptionOfContents() const
    {
	stringstream ss;
	ss << this->Base::descriptionOfContents() << " ";
	ss << "( Name " << _rep_(this->_Name) << " )";
	return ss.str();
    }





CL_LISPIFY_NAME("otherSidePlugName");
CL_DEFMETHOD     core::Symbol_sp Plug_O::otherSidePlugName()
    {
	return DirectionalCoupling_O::otherPlugName(this->_Name);
    }


#ifdef XML_ARCHIVE
    void	Plug_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	if ( node->saving() ) this->getConstitution();
	node->archiveWeakPointer("constitution",this->_WeakConstitution);
#if PRODUCTION_CODE   // FIXME use "name" only and remove the test for "_key"
	node->attribute("name",this->_Name);
#else
	if ( node->loading() )
	{
	    if ( node->hasAttribute("name") )
	    {
		node->attribute("name",this->_Name);
	    } else
	    {
		node->attribute("_key",this->_Name);
	    }
	} else
	{
	    node->attribute("name",this->_Name);
	}
#endif
    }
#endif


    gctools::Vec0<Mate_sp>	_EmptyMateList;

    gctools::Vec0<Mate_sp>::iterator	Plug_O::begin_Mates()
    {
	return _EmptyMateList.end();
    }

    gctools::Vec0<Mate_sp>::iterator	Plug_O::end_Mates()
    {
	return _EmptyMateList.end();
    }














/*
  __BEGIN_DOC(classes.PlugWithMates)
  Defines one or two atoms of this monomer that can be plugged into, a plug name and a
  collection of Mate objects that can act as mates for this plug.
  __END_DOC
*/

/*
  __BEGIN_DOC(classes.PlugWithMates.!class.PlugWithMates)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a PlugWithMates object. PlugWithMatess can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  __END_DOC
*/

string	PlugWithMates_O::descriptionOfContents() const
{
  stringstream ss;
  ss << this->Base::descriptionOfContents() << " ";
  ss << " ( Mates ";
  ss << "#" << this->_Mates.size() << " ";
  gctools::Vec0<Mate_sp>::const_iterator vi;
  for ( vi=this->_Mates.begin(); vi!=this->_Mates.end(); vi++ )
  {
    ss << (*vi)->description() << ", ";
  }
  ss << " )";
  return ss.str();
}


    string PlugWithMates_O::__repr__() const
    {
	stringstream ss;
	ss << "( " << this->className() << " '"<< _rep_(this->getName()) << " :mates ";
	ss << "(list ";
	for ( Mates::const_iterator it=this->_Mates.begin(); it!=this->_Mates.end(); it++ )
	{
          ss << _rep_((*it)) << " ";
	}
	return ss.str();
    }

#ifdef XML_ARCHIVE
    void	PlugWithMates_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attributeIfNotDefault<string>("b0",this->_B0,"");
	node->attributeIfNotDefault<string>("b1",this->_B1,"");
	node->archiveVector0("mates",this->_Mates);
    }
#endif



    bool	PlugWithMates_O::recognizesMateNameOrPdb(core::Symbol_sp name)
    {
        gctools::Vec0<Mate_sp>::iterator	mi;
	for ( mi=this->_Mates.begin(); mi!=this->_Mates.end(); mi++ ) {
	    if ( (*mi)->recognizesNameOrPdb(name) ) return true;
	}
	return false;
    }


/*
  __BEGIN_DOC(classes.OutPlug.!class.OutPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a OutPlug object. OutPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  Outgoing plugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  __END_DOC
*/

    string	OutPlug_O::descriptionOfContents() const
    {
	stringstream ss, shex;
	ss << this->Base::descriptionOfContents() << " ";
	return ss.str();
    }



#ifdef XML_ARCHIVE
    void	OutPlug_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute("stubPivotAtom",this->_StubPivotAtom);
    }
#endif



/*
  __BEGIN_DOC(classes.InPlug.!class.InPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}

  Initialize a InPlug object. InPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  __END_DOC
*/
    string	InPlug_O::descriptionOfContents() const
    {
	return this->Base::descriptionOfContents();
    }


/*
  __BEGIN_DOC(classes.JumpPlug.!class.JumpPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{originFrame:}{Frame::originFrame}

  Initialize a JumpPlug object. JumpPlugs don't make bonds don't have mates but they do have an origin frame that is within the topology.
  __END_DOC
*/

string	JumpPlug_O::descriptionOfContents() const
{
  stringstream ss;
  ss << this->Base::descriptionOfContents() << " ";
  ss << "JumpAtomName[" << this->_JumpAtomName << "]";
  return ss.str();
}


#ifdef XML_ARCHIVE
    void	JumpPlug_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute("jumpAtom",this->_JumpAtomName);
    }
#endif

/*
  __BEGIN_DOC(classes.RingClosingPlug.!class.RingClosingPlug)
  \requiredKeyed{name:}{Text::plugName}
  \requiredKeyed{bond0:}{Text::bond0AtomName}
  \optionalKeyed{bond1:}{Text::bond1AtomName}
  \requiredKeyed{mates:}{Cons::listOfMates}
  \optionalKeyed{exportFrame:}{Frame::exportFrame}
  \requiredKeyed{ringClosingMates:}{Cons::listOfRingClosingMates}

  Initialize a RingClosingPlug object. RingClosingPlugs can have one bond (eg: amide) or two bonds (eg:diketopiperazine).
  RingClosingPlugs export a frame of reference to the next monomer, use \sa{exportFrame} to define this.
  There is an additional list of ringClosingMates that can be attached to this plug without building the mates.
  __END_DOC
*/

#if 0
RingClosingPlug_O::RingClosingPlug_O(const RingClosingPlug_O& p) : RingClosingPlug_O::Base(p)
    {
	LOG(BF("Original %s") % p.description().c_str()  );
        gctools::Vec0<RingClosingMate_sp>::const_iterator vi;
	LOG(BF("Copying %d ring closing mates") % p._Mates.size()  );
	for ( vi=p._RingClosingMates.begin(); vi!=p._RingClosingMates.end(); vi++ ) {
	    LOG(BF("Copied mate") );
	    GC_COPY(RingClosingMate_O, rn, *vi->get() ); // = RP_Copy<RingClosingMate_O>(*vi);
	    this->_RingClosingMates.push_back(rn);
	}
	LOG(BF("Copy constructed %s") % this->description().c_str()  );
    }
#endif


    string	RingClosingPlug_O::descriptionOfContents() const
    {
	stringstream ss, shex;
	ss << this->Base::descriptionOfContents() << " ";
	ss << " ( RingClosingMates ";
	ss << "#" << this->_Mates.size() << " ";
        gctools::Vec0<Mate_sp>::const_iterator vi;
	for ( vi=this->_Mates.begin(); vi!=this->_Mates.end(); vi++ ) {
	    ss << (*vi)->description() << ", ";
	}
	ss << " )";
	return ss.str();
    }


#ifdef XML_ARCHIVE
    void	RingClosingPlug_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->archiveVector0("ringCLosingMates",this->_RingClosingMates);
    }
#endif


















    string Mate_O::__repr__() const
    {
	stringstream ss;
	ss << "( " << this->className();
	{
	    if ( this->getName().notnilp() )
	    {
              ss << " :name " << _rep_(this->getName());
	    }
	}
	ss << _rep_(this->_EntityNames);
	ss << " )";
	return ss.str();

    }

























CL_LISPIFY_NAME("ringClosingMatesAsCons");
CL_DEFMETHOD     core::List_sp RingClosingPlug_O::ringClosingMatesAsCons()
    {_OF();
	return core::Cons_O::createFromVec0(this->_Mates);
    }


    bool RingClosingPlug_O::recognizesRingClosingMate(core::Symbol_sp mateName)
    {_OF();
        gctools::Vec0<Mate_sp>::iterator it;
	for ( it=this->_Mates.begin(); it!=this->_Mates.end(); it++ )
	{
	    if ( (*it)->recognizesMonomerName(mateName)) return true;
	}
	return false;
    }



}; //namespace chem
