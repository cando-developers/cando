       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffStretch.h
 *
 *	Maintains a database of stretch types
 */

#ifndef FFSTR_H
#define	FFSTR_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/adapt/quickDom.h>
#include <cando/adapt/symbolMap.h>
#include <clasp/core/smallMap.h>
#include <cando/chem/ffBaseDb.h>



namespace       chem {



#define	FFStretchWildCardCode	-1

    SMART(FFStretch);
    class FFStretch_O : public FFParameter_O
    {
	LISP_CLASS(chem,ChemPkg,FFStretch_O,"FFStretch",FFParameter_O);
    public:
	void initialize();
        bool fieldsp() const { return true; };
        void fields(core::Record_sp node);
    public:
        core::Symbol_sp	_Type1;
        core::Symbol_sp	_Type2;
	double		_R0_Nanometer;
	double		_Kb_kJPerNanometerSquared;
//	double		_K3;
//	double		_K4;
//	double		_Bci;

    public:
        static FFStretch_sp create_missing(core::Symbol_sp t1, core::Symbol_sp t2) {
          FFStretch_sp stretch = FFStretch_O::create();
          stretch->_Level = unknown;
          stretch->_Type1 = t1;
          stretch->_Type2 = t2;
          return stretch;
        }
    public:
	void setKb_kJPerNanometerSquared(double val);
	void setR0_Nanometer(double val);

	double getR0_Nanometer() const;
	double getR0_Angstrom() const;
	double getKb_kJPerNanometerSquared() const;
	double getKb_kCalPerAngstromSquared() const;
	virtual	string	levelDescription();
	DEFAULT_CTOR_DTOR(FFStretch_O);
    };




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
    inline	string	XmlTag_EstimateStretch() { return "EstimateStretch";};

    class	EstimateStretch {
    public:
      core::Symbol_sp		_ti;
      core::Symbol_sp		_tj;
      double		_rij;
      double		_lnKij;

	adapt::QDomNode_sp	asXml();
	void		parseFromXml(adapt::QDomNode_sp xml);
    };


    inline	string	XmlTag_FFStretchDb() { return "FFStretchDb";};
    class	FFStretchDb_O;


    SMART(FFStretchDb);
    class FFStretchDb_O : public FFBaseDb_O
    {
	LISP_CLASS(chem,ChemPkg,FFStretchDb_O,"FFStretchDb",FFBaseDb_O);

    public:
        bool fieldsp() const { return true; };
        void	fields(core::Record_sp node);
    public:
	typedef	gctools::Vec0<FFStretch_sp>::iterator	FFStretch_spIterator;
        gctools::Vec0<FFStretch_sp>		_Terms;
	adapt::SymbolMap<FFStretch_O>	        _Lookup;
        gctools::SmallMap<core::Symbol_sp,EstimateStretch>	_EstimateStretch; // Not archived

	void	clearEstimateStretch();
	void	addEstimateStretch(core::Symbol_sp ti, core::Symbol_sp tj, double rij, double lnKij );
	void	_addEstimateStretch(const EstimateStretch& es);
	void	add( FFStretch_sp str );
	FFStretch_sp	findTerm(chem::Atom_sp a1, chem::Atom_sp a2 );
//        bool            hasTerm(string t1, string t2);
//        void            cantFind(string t1, string t2);

	FFStretch_spIterator	begin() { return this->_Terms.begin(); };
	FFStretch_spIterator	end() { return this->_Terms.end(); };

	void		initialize();

	DEFAULT_CTOR_DTOR(FFStretchDb_O);
    };



};



TRANSLATE(chem::FFStretch_O);
TRANSLATE(chem::FFStretchDb_O);
#endif

