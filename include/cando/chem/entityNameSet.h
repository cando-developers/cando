/*
    File: entityNameSet.h
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


#ifndef MonomerRecognizer_H
#define MonomerRecognizer_H
#include <clasp/core/common.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
//#include	"conformation.h"
#include <cando/chem/atom.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/adapt/stringList.fwd.h>
#include <clasp/core/symbol.fwd.h>
#include <cando/adapt/symbolList.fwd.h>
#include <cando/chem/entity.h>
#include <cando/chem/statusTracker.h>


#include <cando/chem/chemPackage.h>

#include <cando/adapt/symbolList.fwd.h>// entityNameSet.h wants SymbolList needs symbolList.fwd.h

namespace chem {

SMART(Alias);
SMART(AtomIndexer);
SMART(CandoDatabase);
SMART(SymbolList);

SMART(EntityNameSetBase);
class EntityNameSetBase_O : public Entity_O
{
    LISP_CLASS(chem,ChemPkg,EntityNameSetBase_O,"EntityNameSetBase",Entity_O);
#if INIT_TO_FACTORIES
 public:
    static EntityNameSetBase_sp make(core::List_sp entityNamesCons);
#else
public:
    DECLARE_INIT();
#endif
public:
    void initialize();
    string description() const;
    friend class CandoDatabase_O;
public:
    typedef	adapt::SymbolSet_sp	EntityNames;
protected:
    /*!Ring closing plugs have optional neighbors */
    bool		_Optional;
    adapt::SymbolSet_sp	_EntityNames;
protected:
    void _clear();
    void _removeMonomerName(core::Symbol_sp nm);
public:
    static EntityNameSetBase_sp	create2(core::Symbol_sp nm);
public:
//	virtual void setCandoDatabase(CandoDatabase_sp db);

    virtual adapt::SymbolSet_sp	expandedNameSet();
    virtual RepresentativeList_sp expandedRepresentativeList() const;

    /*! Expand other and merge its names with mine */
    void mergeEntityNames(Entity_sp other);

public: // Lisp creation

public:
CL_LISPIFY_NAME("getMonomerNameWithAtoms");
CL_DEFMETHOD     virtual string getMonomerNameWithAtoms(core::Symbol_sp nm) { return "--- no atom names for Entity_O ---"; };

    void setOptional(bool b) { this->_Optional = b; };
    bool getOptional() { return this->_Optional; };

    void	setGroupNames(adapt::SymbolSet_sp s);


    /*! Return all other monomer names that are not recognized by
     * this EntityNameSetBase
     */
    adapt::SymbolSet_sp	getUnrecognizedMonomerNames();

    /*!Set this EntityNameSetBase up containing one monomer and
     * give it the same name as the monomer
     */
    virtual void	setMonomerNameOrPdb(core::Symbol_sp nm);


    /*! Add the group name and immediatly expand it to monomer names
     */
    virtual void	addMonomerName(core::Symbol_sp nm );
    virtual void	addEntityName(core::Symbol_sp nm ) { this->addMonomerName(nm);};
    virtual void	addMonomerNameNoChecks(core::Symbol_sp nm );


    bool	recognizesMonomerName(core::Symbol_sp name);
    //! Converts nameOrPdb to monomer name and checks if we recognize it
    bool	recognizesNameOrPdb(core::Symbol_sp nameOrPdb );

    adapt::SymbolSet_sp	getUnExpandedNames();

    adapt::SymbolSet_sp	getMonomerNames();
    string		getMonomerNamesAsString();

    //! String together group and monomer names
    virtual core::Symbol_sp getKey();


    /*! Return a StringSet that combines every equivalentName
      from this EntityNameSetBase with every
      equivalentName from the passed EntityNameSetBase.
      This generates a list of monomer pairs that can be
      coupled by a certain coupling.
    */
//    adapt::StringSet_sp combineMonomerNames(string ins, EntityNameSetBase_sp e);



    core::Symbol_sp	getOnlyMonomerName();

    virtual bool isFinalized() { return true; };

    /*! Some subclasses support InterestingAtomAliases
     * This method can be called to figure out which ones do
     */
CL_LISPIFY_NAME("supportsInterestingAtomAliases");
CL_DEFMETHOD     virtual	bool	supportsInterestingAtomAliases() { return false;};
CL_LISPIFY_NAME("setInterestingAtomNamesForMonomerName");
CL_DEFMETHOD     virtual void setInterestingAtomNamesForMonomerName(core::Symbol_sp monomerName, const string& atomIndexerNames) {_OF(); SUBCLASS_MUST_IMPLEMENT();};
CL_LISPIFY_NAME("getInterestingAtomNamesForMonomerName");
CL_DEFMETHOD     virtual string getInterestingAtomNamesForMonomerName(core::Symbol_sp nm) { return "";};
    virtual bool hasInterestingAtomAlias(Alias_sp alias);
    virtual int getInterestingAtomAliasIndex(Alias_sp alias);
    virtual AtomIndexer_sp getAtomIndexerForMonomerName(core::Symbol_sp nm) {_OF();SUBCLASS_MUST_IMPLEMENT();};
    virtual core::List_sp getInterestingAtomAliases();

CL_LISPIFY_NAME("getInterestingAtomAliasesAsString");
CL_DEFMETHOD     virtual string getInterestingAtomAliasesAsString() { return "";};

    adapt::StringList_sp	getMonomerNamesOrdered();
    adapt::StringList_sp	getUnrecognizedMonomerNamesOrdered();



    /*! Expand my EntityNames to terminal EntityNames */
    void expandToTerminalEntityNames();


    /*! Contract the Entity names */
    void contractEntityNames(adapt::SymbolSet_sp entityNameSets);


CL_LISPIFY_NAME("testEntityNameSetBase");
CL_DEFMETHOD     string	testEntityNameSetBase() { return "testEntityNameSetBase"; };

    EntityNameSetBase_O( const EntityNameSetBase_O& emr );

    DEFAULT_CTOR_DTOR(EntityNameSetBase_O);
};




/*! All EntityNameSet adds to EntityNameSetBase is that it stores a name for the EntityNameSet
 */
SMART(EntityNameSet);
class EntityNameSet_O : public EntityNameSetBase_O
{
    friend class CandoDatabase_O;
    LISP_CLASS(chem,ChemPkg,EntityNameSet_O,"EntityNameSet",EntityNameSetBase_O);
#if INIT_TO_FACTORIES
 public:
    static EntityNameSet_sp make(core::Symbol_sp name);
#else
    DECLARE_INIT();
#endif
public:
	void initialize();
	string description() const;
public:
	static EntityNameSet_sp	create2(core::Symbol_sp nm);
public:
    core::Symbol_sp		_Name;
public:
CL_LISPIFY_NAME("setName");
CL_DEFMETHOD     void setName(core::Symbol_sp nm) { this->_Name = nm;};
CL_LISPIFY_NAME("getName");
CL_DEFMETHOD 	core::Symbol_sp getName() const { return this->_Name;};
	EntityNameSet_O( const EntityNameSet_O& emr );

	DEFAULT_CTOR_DTOR(EntityNameSet_O);
};



SMART(EntityNameSetWithCap);
class EntityNameSetWithCap_O : public EntityNameSetBase_O
{
    LISP_CLASS(chem,ChemPkg,EntityNameSetWithCap_O,"EntityNameSetWithCap",EntityNameSetBase_O);
public:
    void	initialize();
    string description() const;
//	virtual void setCandoDatabase(CandoDatabase_sp db);
    /*! Add the group name and immediately expand it to
     * a list of monomers */
    void	addGroupName(core::Symbol_sp nm );
    bool	recognizesGroupName(core::Symbol_sp nm);


    EntityNameSetWithCap_O( const EntityNameSetWithCap_O& emr );
    DEFAULT_CTOR_DTOR(EntityNameSetWithCap_O);
};




};










TRANSLATE(chem::EntityNameSetBase_O);
TRANSLATE(chem::EntityNameSet_O);
TRANSLATE(chem::EntityNameSetWithCap_O);
#endif
