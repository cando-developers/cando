/*
    File: aggregate.h
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





#ifndef	AGGREGATE_H
#define	AGGREGATE_H

#include <clasp/core/common.h>
#include <string>
#include <cando/chem/matter.h>
#include <cando/chem/atomIdMap.fwd.h>
#include <cando/chem/molecule.h>



#include <cando/chem/chemPackage.h>

namespace chem {

class AtomInfo {
public:
  Atom_sp         _Atom;
  Residue_sp      _Residue;
  core::Symbol_sp _Type;
  size_t          _MoleculeIndex;
  size_t          _ResidueIndex;
};

class ResidueOut {
public:
  Molecule_sp    mol;
  Residue_sp     res;
  uint           firstAtom;
};

SMART(Restraint);
SMART(Lisp);

FORWARD(BoundingBox);
class  BoundingBox_O : public  core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,BoundingBox_O,"BoundingBox", core::CxxObject_O );
//    DECLARE_SERIALIZE();
public:
  Vector3 _Center;
  Vector3 _Widths;
  Vector3 _AnglesDegrees;
  // for optimization
  double  _x_rwidth;
  double  _y_rwidth;
  double  _z_rwidth;
public:
  bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
public:
  BoundingBox_O() {};
  static BoundingBox_sp make(core::List_sp widths, core::T_sp angles_degrees, core::T_sp center);
public:
  Vector3 get_bounding_box_widths() const;
  Vector3 get_bounding_box_rwidths() const;
  Vector3 get_bounding_box_angles_degrees() const;
  Vector3 get_bounding_box_center() const;
  bool cuboidp() const;
  Vector3 min_corner() const;
  Vector3 max_corner() const;
  Vector3 normalize_position(const Vector3& pos) const;
  double get_x_width() const;
  double get_y_width() const;
  double get_z_width() const;
  double get_x_rwidth() const { return this->_x_rwidth; }
  double get_y_rwidth() const { return this->_y_rwidth; }
  double get_z_rwidth() const { return this->_z_rwidth; }
  double get_x_angle_degrees() const;
  double get_y_angle_degrees() const;
  double get_z_angle_degrees() const;
  core::T_mv get_cuboid_rsize() const;
  std::string __repr__() const;


  double distance_squared_between_two_points(const Vector3& v1, const Vector3& v2);
  double distance_squared_between_two_atoms(Atom_sp a1, Atom_sp a2);

  bool intersects(BoundingBox_sp other);

  BoundingBox_O(const Vector3& widths, const Vector3& angles_degrees, const Vector3& center);
};

class  Aggregate_O : public  Matter_O
{
    LISP_CLASS(chem,ChemPkg,Aggregate_O,"Aggregate", Matter_O );
//    DECLARE_SERIALIZE();
public:
  core::T_sp _ForceFieldName;
public:
    void initialize();
    friend	class	Loop;
  BoundingBox_sp _BoundingBox;
public:
  bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
public:
  static Aggregate_sp make(MatterName name, core::List_sp molecules = nil<core::T_O>());

public:
    static Aggregate_sp	open(const string& name);
public:

    /*! Build an AtomIdMap that maps moleculeId/residueId/atomId to atoms
     Call this after you make any changes to the structure of the Aggregate.
    This is automatically called whenever an Aggregate is loaded. */
    virtual AtomIdMap_sp buildAtomIdMap() const;


    virtual Atom_sp atomWithAtomId(const AtomId& atomId) const;

  BoundingBox_sp boundingBox() const;
  bool boundingBoxBoundP() const;
  void setBoundingBox(BoundingBox_sp bounding_box);
  void makUnboundBoundingBox();

    /*! Lookup an Atom using the atomId */
//    Atom_sp lookupAtom(const AtomId& atomId) const;


    char	getMatterType() { return AGGREGATE_CLASS; };

  contentIterator begin_molecules() { return this->_Contents.begin(); };
  contentIterator end_molecules() {return this->_Contents.end(); };

    Atom_sp		atomWithId(int i);

//    virtual bool equal(core::T_sp obj) const;
    virtual void	transferCoordinates(Matter_sp other);

    virtual Matter_mv addMatter(Matter_sp matter);
CL_LISPIFY_NAME("addMolecule");
CL_DEFMETHOD     void	addMolecule( Molecule_sp a ) {this->addMatter(a);};
    Aggregate_sp	addMoleculeRetainId( Molecule_sp a ) {this->addMatterRetainId(Matter_sp(a)); return this->sharedThis<Aggregate_O>();};
    Molecule_sp      firstMolecule();
    MatterName     firstMoleculeName();

//    Aggregate_sp	shallowCopy();
//    Aggregate_sp	deepCopy() const;

    core::List_sp atomsWithChimeraSpecifications(const string& name);
    Atom_sp		firstAtomWithName( MatterName name );

    void		perturbAtomPositions(double dist);

    void		removeMolecule( Molecule_sp mol );
//	void		dump();

    bool		testForVdwOverlap(Aggregate_sp other);

    virtual uint	numberOfAtoms( );
    virtual bool isAggregate() { return true;};

    /*! Return a Cons with all molecules */
CL_LISPIFY_NAME("molecules");
CL_DEFMETHOD     core::List_sp molecules() { return this->contentsAsList(); };

    /*! When a structure is loaded and all the atoms are dumped into the same
     * molecule and same residue then we may want to separate out the molecules
     * defined by their connectivity.
     *
     * This function identifies the separate molecules using spanning trees
     * and creates a separate molecule for each of them.
     *
     * Return the number of molecules found.
     */
    uint separateMolecules();

    /*! Sometimes we want to identify a molecule that contains an atom with
     * a specific name.
     */
    Molecule_sp firstMoleculeWithAtomNamed(MatterName name);


    void		writeToFile(const string& fileName);
CL_LISPIFY_NAME("saveAs");
CL_DEFMETHOD     void		saveAs(const string& fileName) { _G(); this->writeToFile(fileName);};


    bool	canRender() { return true;}
#ifdef RENDER
    geom::Render_sp	rendered(core::List_sp opt);
#endif


    virtual string	subMatter() { return Molecule_O::static_className(); };
    virtual	string	description() const { stringstream ss; ss << "aggregate("<<_rep_(this->getName())<<")@"<<std::hex<<this<<std::dec; return ss.str();};


    void		duplicate(const Aggregate_O * x );


public:
    virtual Matter_sp	copy(core::T_sp new_to_old);

    virtual Matter_sp copyDontRedirectAtoms(core::T_sp new_to_old);
	virtual void redirectAtoms();

  core::T_sp force_field_name() const;
  void setf_force_field_name(core::T_sp name);


  Aggregate_O() : _ForceFieldName(kw::_sym_default), _BoundingBox(unbound<BoundingBox_O>()) {};
};


};
TRANSLATE(chem:: Aggregate_O);
#endif
