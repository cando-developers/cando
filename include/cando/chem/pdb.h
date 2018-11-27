/*
    File: pdb.h
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
  FORWARD(Atom);
    struct	AtomPdbRec
    {
	string _line;
	// pdb stuff
	string	_recordName;
	int	_serial;
	core::Symbol_sp _name;
	string	_altLoc;
	core::Symbol_sp _resName;
	string _chainId;
	int	_resSeq;
	string	_iCode;
	double	_x;
	double	_y;
	double	_z;
	double	_occupancy;
	double	_tempFactor;
	string	_element;
	string	_charge;
	// my stuff
	Atom_sp	_atom;
	int	_residueIdx;
	int	_moleculeIdx;

        AtomPdbRec() : _atom(_Nil<Atom_O>()) {};
	void write(core::T_sp stream);
	virtual ~AtomPdbRec() {};
	void parse(const string& line);
	Atom_sp createAtom();
    };


};
namespace chem
{

    SMART(Aggregate);
    SMART(Matter);

#if 0
    __ BEGIN_CLASS_DEFINITION(ChemPkg,PdbMonomerConnectivity_O,PdbMonomerConnectivity,core::T_O)
	private:
	string _PdbName;
    geom::ObjectList_sp _HetNames;
gc::SmallMap<core::Symbol_sp,adapt::StringSet_sp> _Connections;
public:

    __END_CLASS_DEFINITION(PdbMonomerConnectivity_O)


    __ BEGIN_CLASS_DEFINITION(ChemPkg,PdbMonomerDatabase_O,PdbMonomerDatabase,core::T_O)
    __END_CLASS_DEFINITION(PdbMonomerDatabase_O)
#endif





    SMART(PdbReader );
    class PdbReader_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,PdbReader_O,"PdbReader",core::CxxObject_O);
    public:
	void	initialize();
    private:
    public:
	static Aggregate_sp loadPdb(core::T_sp fileName);
	static Aggregate_sp loadPdbConnectAtoms(core::T_sp fileName);
        static Aggregate_sp loadPdbFromStreamConnectAtoms(core::T_sp stream);
    public:

	Aggregate_sp	parse(core::T_sp stream);

	PdbReader_O( const PdbReader_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(PdbReader_O);
    };




    SMART(PdbWriter );
    class PdbWriter_O : public core::CxxObject_O
    {
      LISP_CLASS(chem,ChemPkg,PdbWriter_O,"PdbWriter",core::CxxObject_O);
#if INIT_TO_FACTORIES
    public:
      static PdbWriter_sp make(core::T_sp fileName);
#else
      DECLARE_INIT();
#endif
    public:
      void	initialize();
    private:
      core::T_sp      _Out;
    public:
      static void savePdb(Matter_sp matter, core::T_sp fileName);
    public:

      void open(core::T_sp fileName);
      void close();

      void write(Matter_sp matter);
      void writeModel(Matter_sp matter, int model);

      PdbWriter_O( const PdbWriter_O& ss ); //!< Copy constructor

    PdbWriter_O() :  _Out(_Nil<core::T_O>()) {};
      virtual ~PdbWriter_O() {};
    };

};
TRANSLATE(chem::PdbReader_O);
TRANSLATE(chem::PdbWriter_O);
#endif //]
