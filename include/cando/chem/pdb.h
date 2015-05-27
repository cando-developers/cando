#ifndef	PdbReader_H //[
#define PdbReader_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/chemPackage.h>
#include <cando/adapt/stringSet.fwd.h>


namespace chem
{

    SMART(Aggregate);
    SMART(Matter);

#if 0
    __ BEGIN_CLASS_DEFINITION(ChemPkg,PdbMonomerConnectivity_O,PdbMonomerConnectivity,core::T_O)
	private:
	string _PdbName;
    geom::ObjectList_sp _HetNames;
    adapt::SymbolMap<adapt::StringSet_O> _Connections;
public:

    __END_CLASS_DEFINITION(PdbMonomerConnectivity_O)


    __ BEGIN_CLASS_DEFINITION(ChemPkg,PdbMonomerDatabase_O,PdbMonomerDatabase,core::T_O)
    __END_CLASS_DEFINITION(PdbMonomerDatabase_O)
#endif





    SMART(PdbReader );
    class PdbReader_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,PdbReader_O,"PdbReader");
    public:
//	void	archive(core::ArchiveP node);
	void	initialize();
    private:
    public:
	static Aggregate_sp loadPdb(const string& fileName);
	static Aggregate_sp loadPdbConnectAtoms(const string& fileName);
    public:

	Aggregate_sp	parse(const string& fileName);

	PdbReader_O( const PdbReader_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(PdbReader_O);
    };




    SMART(PdbWriter );
    class PdbWriter_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,PdbWriter_O,"PdbWriter");
#if INIT_TO_FACTORIES
    public:
	static PdbWriter_sp make(string fileName);
#else
	DECLARE_INIT();
#endif
    public:
	void	initialize();
    private:
        core::T_sp      _Out;
    public:
	static void savePdb(Matter_sp matter, const string& fileName);
    public:

	void open(const string& fileName);
	void close();

	void write(Matter_sp matter);
	void writeModel(Matter_sp matter, int model);

	PdbWriter_O( const PdbWriter_O& ss ); //!< Copy constructor

	PdbWriter_O() :  _Out(_Nil<core::T_O>()) {};
	virtual ~PdbWriter_O();
    };


};
TRANSLATE(chem::PdbReader_O);
TRANSLATE(chem::PdbWriter_O);
#endif //]
