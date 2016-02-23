/*
    File: rotamer.cc
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
#include <cando/chem/rotamer.h>
#include <clasp/core/wrappers.h>




namespace chem
{



#ifdef XML_ARCHIVE
void	RotamerAtom::archive(core::ArchiveP node)
{
    node->attribute("d",this->_Distance);
    node->attribute("a",this->_Angle);
    node->attribute("t",this->_Torsion);
    LOG(BF("Serialized RotamerAtom@%X = %lf, %lf, %lf") % this % this->_Distance % this->_Angle % this->_Torsion );
};

#endif







// ----------------------------------------------------------------------
//







#ifdef XML_ARCHIVE
    void Rotamer_O::archiveBase(core::ArchiveP node)
    {
	if ( node->saving() ) {
	    uint sz = this->_RotamerAtoms.size();
	    node->attribute("size",sz);
	    int i;
	    core::ArchiveP vnode;
	    RotamerAtoms::iterator vi;
	    for ( i=0,vi=this->_RotamerAtoms.begin(); vi!=this->_RotamerAtoms.end(); i++,vi++ ) {
		vnode = node->createChildNode("V");
		vnode->attribute("i",i);
		(*vi).archive(vnode);
	    }
	} else {
	    uint sz=0, i=0;
	    node->attribute("size",sz);
	    this->_RotamerAtoms.resize(sz);
	    core::VectorArchiveNodes::iterator vi;
	    i = 0;
	    for ( vi=node->begin_Children(); vi!=node->end_Children(); vi++ ) {
		(*vi)->setRecognized(true);
		(*vi)->attribute("i",i);
		this->_RotamerAtoms[i].archive(*vi);
		i++;
	    }
	}
    }
#endif


    void Rotamer_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    

}; /* chem */
