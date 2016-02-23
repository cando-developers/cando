/*
    File: atomIdMap.cc
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
#include <clasp/core/environment.h>
#include <cando/chem/atomId.h>
#include <cando/chem/atomIdMap.h>
#include <cando/chem/atom.h>
#include <clasp/core/wrappers.h>
namespace chem
{

// ----------------------------------------------------------------------
//




CL_LISPIFY_NAME("lookup-atom");
CL_DEFMETHOD     Atom_sp AtomIdToAtomMap_O::lookupAtom(AtomId_sp atomId) const
    {_OF();
	return this->_AtomIdMap[atomId->_AtomId];
    }


#ifdef XML_ARCHIVE
    void AtomIdToAtomMap_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	IMPLEMENT_ME();
    }
#endif


    void AtomIdToAtomMap_O::initialize()
    {_OF();
        this->Base::initialize();
    }


    void AtomIdToAtomMap_O::resize(int numMols)
    {_OF();
	this->_AtomIdMap.resize(numMols);
    }

    void AtomIdToAtomMap_O::resize(int mol, int numRes)
    {_OF();
	this->_AtomIdMap.resize(mol,numRes);
    }

    void AtomIdToAtomMap_O::resize(int mol, int res, int numAtoms)
    {_OF();
	this->_AtomIdMap.resize(mol,res,numAtoms);
    }
    void AtomIdToAtomMap_O::set(AtomId const& atomId, Atom_sp atom)
    {_OF();
	this->_AtomIdMap[atomId] = atom;
    }

    

}; /* chem */
