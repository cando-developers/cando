/*
    File: stereochemistry.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//



/*
 *	stereochemistry.h
 *
 *
 */

#ifndef Stereochemistry_H
#define Stereochemistry_H
#include <clasp/core/common.h>

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/adapt/symbolSet.fwd.h>
#include <cando/geom/vector2.h>
#include <cando/geom/vector3.h>
#include <cando/chem/entity.h>


#include <cando/chem/chemPackage.h>


namespace chem
{
  SMART(Atom);
  SMART(Residue);
  SMART(Topology);
  SMART(ComplexRestraint);

  SMART(StereoConfiguration);
  class StereoConfiguration_O : public core::CxxObject_O
  {
    CL_DOCSTRING(R"(Store the stereochemical configuration :R or :S for a single atom represented by the atoms name)");
    LISP_CLASS(chem,ChemPkg,StereoConfiguration_O,"StereoConfiguration",core::CxxObject_O);
  public:
    static StereoConfiguration_sp make(core::Symbol_sp atomName, core::Symbol_sp config);
  private:
    MatterName	_AtomName;
    core::Symbol_sp 	_Configuration;
  public:
    StereoConfiguration_O() {};
  StereoConfiguration_O(MatterName atomName, core::Symbol_sp configuration) : _AtomName(atomName), _Configuration(configuration) {};
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    CL_LISPIFY_NAME(stereochemical_permuations);
    CL_DEF_CLASS_METHOD static	core::List_sp stereochemical_permutations(uint numberOfCenters);
    CL_LISPIFY_NAME(stereo_configuration_create_multiple);
    CL_DEF_CLASS_METHOD static core::List_sp create_multiple(core::List_sp atomNames, core::List_sp configurations);
  public:
    CL_DEFMETHOD MatterName getAtomName()	{return this->_AtomName;};
    CL_DEFMETHOD void setAtomName(MatterName n) {this->_AtomName = n;};
    CL_DEFMETHOD core::Symbol_sp getConfiguration()	{return this->_Configuration;};
    CL_DEFMETHOD void setConfiguration(core::Symbol_sp p) {this->_Configuration = p;};
    int		getMoeConfiguration();
  };

/*! Save a monomer name and pdb name pair
 */
  SMART(Stereoisomer);
  class Stereoisomer_O : public Entity_O {
    CL_DOCSTRING(R"(Stores the stereochemical configuration of a collection of atoms and a name for the stereoisomer)");
    LISP_CLASS(chem,ChemPkg,Stereoisomer_O,"Stereoisomer",Entity_O);
  public:
    static Stereoisomer_sp make(core::Symbol_sp name, core::Symbol_sp pdb, core::Integer_sp stereoisomer_index, core::List_sp configs);
  public:
    core::Symbol_sp	_Name;
    core::Symbol_sp	_Pdb;
    core::Symbol_sp	_Enantiomer;
    core::Integer_sp    _StereoisomerIndex;
	    // Add configurations here
    gctools::Vec0<StereoConfiguration_sp>	_Configurations;
  public:
    Stereoisomer_O() {};
  Stereoisomer_O(core::Symbol_sp name, core::Symbol_sp pdbname, core::Integer_sp stereoisomer_index) : _Name(name), _Pdb(pdbname), _Enantiomer(name), _StereoisomerIndex(stereoisomer_index) {};
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:		// local implementations of CandoDatabaseDependent functions
    bool	isTerminalName() { return true; };
    adapt::SymbolSet_sp	expandedNameSet();
    RepresentativeList_sp expandedRepresentativeList() const;
    core::Symbol_sp getConfigurationForCenter(core::Symbol_sp centerName );
  public:
    CL_DEFMETHOD core::Symbol_sp getName() const {return this->_Name;};
    CL_DEFMETHOD void setName(core::Symbol_sp n) {this->_Name = n;};
    CL_DEFMETHOD core::Integer_sp getStereoisomerIndex() const {return this->_StereoisomerIndex;};
    CL_DEFMETHOD void setStereoisomerIndex(core::Integer_sp n) {this->_StereoisomerIndex = n;};
    CL_DEFMETHOD core::Symbol_sp getPdb() {return this->_Pdb;};
    CL_DEFMETHOD void setPdb(core::Symbol_sp p) {this->_Pdb = p;};
    core::Symbol_sp	getEnantiomer()	{return this->_Enantiomer;};
    void	setEnantiomer(core::Symbol_sp p) {this->_Enantiomer = p;};
    void	addStereoConfiguration(StereoConfiguration_sp sc) {this->_Configurations.push_back(sc);};
    string __repr__() const;
    gctools::Vec0<StereoConfiguration_sp>::iterator _Configurations_begin() { return this->_Configurations.begin();};
    gctools::Vec0<StereoConfiguration_sp>::iterator _Configurations_end() { return this->_Configurations.end();};
  };

  SMART(StereoInformation);
  class StereoInformation_O : public core::CxxObject_O
  {
    CL_DOCSTRING(R"(Store a list of stereoisomers and a mapping of names to stereoisomers. The names are compared using EQ and they can be symbols or fixnums that represent absolute stereochemistry using a canonical order of stereocenters that si defined by the individual stereoisomers and must be consistent across them.)");
    LISP_CLASS(chem,ChemPkg,StereoInformation_O,"StereoInformation",core::CxxObject_O);
  public:
    CL_LISPIFY_NAME(make_stereoinformation);
    CL_DEF_CLASS_METHOD static StereoInformation_sp make(core::List_sp stereoisomers, core::List_sp restraints);
  private:
    gctools::Vec0<Stereoisomer_sp>		_Stereoisomers;
    gc::SmallMap<core::T_sp,Stereoisomer_sp>	_NameOrPdbToStereoisomer;
    gctools::Vec0<ComplexRestraint_sp>	_ComplexRestraints;
  public:
    StereoInformation_O() {};
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    gctools::Vec0<Stereoisomer_sp>::iterator begin_Stereoisomers() { return this->_Stereoisomers.begin();};
    gctools::Vec0<Stereoisomer_sp>::iterator end_Stereoisomers() { return this->_Stereoisomers.end();};

    gctools::Vec0<Stereoisomer_sp>::const_iterator begin_Stereoisomers() const { return this->_Stereoisomers.begin();};
    gctools::Vec0<Stereoisomer_sp>::const_iterator end_Stereoisomers() const { return this->_Stereoisomers.end();};


#if 0
    gctools::Vec0<O_ProChiralCenter>::iterator begin_ProChiralCenters() { return this->_ProChiralCenters.begin();};
    gctools::Vec0<O_ProChiralCenter>::iterator end_ProChiralCenters() { return this->_ProChiralCenters.end();};
#endif
    gctools::Vec0<ComplexRestraint_sp>::iterator begin_ComplexRestraints() { return this->_ComplexRestraints.begin();};
    gctools::Vec0<ComplexRestraint_sp>::iterator end_ComplexRestraints() { return this->_ComplexRestraints.end();};

    CL_DEFMETHOD core::List_sp stereoisomersAsList() { return core::Cons_O::createFromVec0(this->_Stereoisomers); };

    void validate();

    void addStereoisomer(Stereoisomer_sp s);
//    void addProChiralCenter(RPProChiralCenter s);
    Stereoisomer_sp	getStereoisomer(core::T_sp nameOrPdb);

	    //! Return true if this StereoInformation recognizes the name or pdb name
    bool	recognizesNameOrPdb(core::T_sp nm)
    {
      return ( this->_NameOrPdbToStereoisomer.contains(nm) != 0);
    };
	    //! The name of the monomer with the nameOrPdb
    core::T_sp	nameFromNameOrPdb(core::T_sp np)
    {
      Stereoisomer_sp i;
      try {i = this->_NameOrPdbToStereoisomer.get(np);}
      catch (...) {CELL_ERROR(np);};
      ASSERTNOTNULL(i);
      return i->getName();
    };
	    //! The pdb_name of the monomer with the nameOrPdb
    core::T_sp	pdbFromNameOrPdb(core::T_sp np)
    {
      Stereoisomer_sp i;
      try {i = this->_NameOrPdbToStereoisomer.get(np);}
      catch (...) {CELL_ERROR(np);};
      ASSERTNOTNULL(i);
      return i->getPdb();
    };

	    //! Return the names of Monomers that can be created from this constitution
    adapt::StringList_sp getMonomerNamesAsStringList();
    adapt::SymbolSet_sp getMonomerNamesAsSymbolSet();
	    //! Return the pdb names of Monomers that can be created from this constitution
    adapt::StringList_sp getPdbNamesAsStringList();
  };


};
#endif
