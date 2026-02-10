// mmcif_writer.h  and mmcif_writer.cc  (single file for brevity)
//
// A very thin mmCIF writer that mirrors the API of PdbWriter_O.
//  • open(pathDesignator)   – open a CLASP stream for output
//  • close()                – write trailing '#' and close the stream
//  • write(matter, metadata) – write a complete mmCIF block
//  • writeModel(matter, model, metadata) – wrap write() with a data_### block per model
//
// NOTE:
//   • Only the essential _cell, _symmetry and _atom_site categories are emitted.
//   • Connectivity (_struct_conn) and secondary‑structure tables are *not* generated yet.
//   • Uses the same helper code as PdbWriter_O (Loop, AtomPdbRec, etc.) so the header
//     list is almost identical.  This should compile & link anywhere PdbWriter_O does.
//
//   Christian Schafmeister, 2025‑07‑19
//
#ifndef MmcifWriter_H
#define MmcifWriter_H

#include <iomanip>
#include <cctype>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>

#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/pathname.h>
#include <clasp/core/hashTableEqual.h>
#include <clasp/core/string.h>
#include <clasp/core/lispStream.h>

#include <cando/chem/matter.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>

namespace chem
{

  SMART(MmcifWriter);

  /*!
   * A minimal writer for PDBx/mmCIF files.
   * The interface mirrors PdbWriter_O so existing callers can switch by changing
   * the constructor.
   */
  class MmcifWriter_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,MmcifWriter_O,"MmcifWriter",core::CxxObject_O);
  public:
    static core::HashTableEqual_sp saveMmcif(Matter_sp matter, core::T_sp fileName, core::List_sp metadata );
    /*! factory */
    static MmcifWriter_sp make(core::T_sp pathDesignator)
    {
      auto me = gctools::GC<MmcifWriter_O>::allocate();
      me->open(pathDesignator);
      return me;
    }

    /*! Create a fresh CLASP stream for output */
    void open(core::T_sp pathDesignator);

    /*! Emit final END record and close stream */
    void close();

    core::HashTableEq_sp collectResidues(Matter_sp matter );
    void writeResidue(Residue_sp residue);
    void writeAllResidues(core::HashTableEq_sp residues);

    /*!
     *  Write one mmCIF data block.
     *  @param matter   Aggregate/Molecule/Residue hierarchy to serialize
     *  @param metadata (&key :cryst1 (a b c alpha beta gamma space-group Z))
     *
     *  Returns a hash‑table mapping long residue names -> three‑letter comp IDs
     *  exactly like PdbWriter_O does, so downstream code that relies on the
     *  mapping will continue to work.
     */
    core::HashTableEqual_sp write(Matter_sp matter,
                                  core::List_sp metadata=nil<core::T_O>());

    /*! convenience wrapper that prefixes the block with "data_<model>"  */
    void writeModel(Matter_sp matter, int model,
                    core::List_sp metadata=nil<core::T_O>());

  MmcifWriter_O() : _Out(nil<core::T_O>()) {}
    virtual ~MmcifWriter_O() {}

  private:
    core::T_sp _Out;

    struct AtomCifRec {
      core::Symbol_sp _name;
      std::string     _element;
      core::Symbol_sp _resName;
      std::string     _pdbResName;
      char            _chainId;
      int             _resSeq;
      double          _x,_y,_z;
      double          _occupancy;
      double          _biso;
      Atom_sp         _atom;

    AtomCifRec() : _atom(nil<Atom_O>()), _chainId('A'),
        _resSeq(0),_x(0),_y(0),_z(0),
        _occupancy(1.0), _biso(1.0) {}
    };

    /*! helper that collects atoms, assigns serial numbers, builds residue map */
    size_t collectAtoms(Matter_sp matter,
                        std::vector<AtomCifRec>& atoms,
                        std::map<std::string,std::string>& resNameMap,
                        char chainIdStart='A');

  };

}; // namespace chem

#endif /* MmcifWriter_H */

