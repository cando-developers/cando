

#include <clasp/core/common.h>
//#include "core/archive.h"
//#include "core/archiveNode.h"
#include <cando/chem/atom.h>
#include <cando/chem/indirectAtomCoordinateReference.h>

namespace chem
{

    IndirectAtomCoordinateReference::IndirectAtomCoordinateReference()
	: _SharedAtom(_Nil<Atom_O>()), _IndexTimes3(UndefinedUnsignedInt) {};

#ifdef XML_ARCHIVE
    void IndirectAtomCoordinateReference::archive(core::ArchiveP node)
    {
	node->attribute("atom",this->_SharedAtom);
	node->attribute("idx3",this->_IndexTimes3);
    };
#endif
};
