#ifndef	VirtualAtom_H //[
#define VirtualAtom_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/atom.h>
#include <clasp/core/executables.h>
#include <cando/chem/chemPackage.h>


namespace chem {


SMART(CalculatePosition);

SMART(VirtualAtom );
class VirtualAtom_O : public Atom_O
{
    LISP_BASE1(Atom_O);
    LISP_CLASS(chem,ChemPkg,VirtualAtom_O,"VirtualAtom");
public:
//    void	archive(core::ArchiveP node);
	void	initialize();
private:
		//! Calculate position using a procedure
	CalculatePosition_sp	_CalculatePosition;
public:
    static VirtualAtom_sp create(MatterName name, CalculatePosition_sp proc);
public:

	void calculatePosition();

	virtual bool	equal(core::T_sp obj) const;

	virtual bool isVirtual() const { return true;};


			// Constructor stuff
	virtual Matter_sp copy();
//	virtual Atom_sp	copyDropExternalResidueBonds();


	VirtualAtom_O( const VirtualAtom_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(VirtualAtom_O);
};


};

TRANSLATE(chem::VirtualAtom_O);

#endif //]
