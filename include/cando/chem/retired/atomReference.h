#define	DEBUG_LEVEL_FULL     
       
#ifndef	AtomReference_H //[
#define AtomReference_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/candoBase/vector3.h>
#include <clasp/core/executables.h>


#include <cando/chem/chemPackage.h>


namespace chem 
{


    SMART(PointProvider );
    SMART(ScorerState);

    SMART(AtomReferenceBase);
    class AtomReferenceBase_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,AtomReferenceBase_O,"AtomReferenceBase");
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
	AtomReferenceBase_O(const core::Lisp_sp& lisp);
	AtomReferenceBase_O( const AtomReferenceBase_O& ss ); //!< Copy constructor
	DEFAULT_CTOR_DTOR(AtomReferenceBase_O);
    };



    SMART(AtomReference);
    class AtomReference_O : public AtomReferenceBase_O
    {
	LISP_BASE1(AtomReferenceBase_O);
	LISP_CLASS(chem,ChemPkg,AtomReference_O,"AtomReference");
#if INIT_TO_FACTORIES
    public:
	static AtomReference_sp make(const string& alias);
#else
	DECLARE_INIT();	
#endif
//	void archiveBase(core::ArchiveP node);
    private:
	string	_AtomAlias;
    public:
	bool equal(core::T_sp obj) const;
	void	initialize();
	DEFAULT_CTOR_DTOR(AtomReference_O);
    };



};
TRANSLATE(chem::AtomReferenceBase_O);
TRANSLATE(chem::AtomReference_O);
#endif //]
