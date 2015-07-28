       
       
#ifndef	Alias_H //[
#define Alias_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <clasp/core/executables.h>


#include <cando/chem/chemPackage.h>


namespace chem
{


    SMART(ScorerState);

    SMART(Alias);
    class Alias_O : public core::CxxObject_O
    {
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,Alias_O,"Alias");
#if INIT_TO_FACTORIES
    public:
	static Alias_sp make(core::Symbol_sp monomerSym, core::Symbol_sp atomSym);
#else
	DECLARE_INIT();
#endif
    public:
    public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
    public:
	static Alias_sp create(core::Symbol_sp monomerAlias, core::Symbol_sp atomAlias )
	{
	    GC_ALLOCATE(Alias_O, a ); // RP_Create<Alias_O>(e);
	    a->setMonomerAlias(monomerAlias);
	    a->setAtomAlias(atomAlias);
	    return a;
	}
	static core::T_sp prim_alias(core::Function_sp, core::List_sp args);
    private:
	core::Symbol_sp	_MonomerAlias;
	core::Symbol_sp	_AtomAlias;
    private: // Do not archive indices
	uint	_MonomerId;
	uint	_AtomId;
    public:

	core::Symbol_sp getMonomerAlias() { return this->_MonomerAlias;};
	void	setMonomerAlias(core::Symbol_sp a) { this->_MonomerAlias = a;};
	core::Symbol_sp getAtomAlias() { return this->_AtomAlias;};
	void	setAtomAlias(core::Symbol_sp a) { this->_AtomAlias = a;};


#if 0
	/*! If the alchemist oligomer changes then update this
	 * alias to the correct monomer index
	 */
	void oligomerChanged(PointProvider_sp alchemist);
	/*! If the alchemist oligomer changes then update this
	 * alias to the correct atom index
	 */
	void sequenceChanged(PointProvider_sp alchemist);
	/*! Get the atom position for this alias
	 */
//	Vector3 getAtomPosition(PointProvider_sp alchemist);

	/*! Get the transformed atom position for this alias
	 */
	Vector3& getUntransformedAtomPosition(PointProvider_sp alchemist);
#endif


	virtual bool equal(core::T_sp obj) const;

	string description() const;
	string __repr__() const { return this->_MonomerAlias->__repr__()+"@"+this->_AtomAlias->__repr__();};
	Alias_O(const core::Lisp_sp& lisp);
	Alias_O( const Alias_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(Alias_O);
    };




};
TRANSLATE(chem::Alias_O);
#endif //]
