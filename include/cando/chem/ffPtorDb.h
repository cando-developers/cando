       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffPtor.h
 *
 *	Maintains a database of ptor terms
 */

#ifndef FFPTOR_H
#define	FFPTOR_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/symbolMap.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/ffBaseDb.h>
#include <clasp/core/smallMap.h>


#include <cando/chem/chemPackage.h>


namespace chem
{


#define	FFPtorWildCardCode	-1


    SMART(FFPtor);
    class FFPtor_O : public FFParameter_O
    {
	LISP_BASE1(FFParameter_O);
	LISP_CLASS(chem,ChemPkg,FFPtor_O,"FFPtor");
    public:
	static int const MaxPeriodicity = 6;
    public:
	void initialize();
#ifdef XML_ARCHIVE
	void	archiveBase(core::ArchiveP node);
#endif

    public:
	core::Symbol_sp	_T1;
	core::Symbol_sp _T2;
	core::Symbol_sp _T3;
	core::Symbol_sp _T4;
	bool		_hasPeriodicity[MaxPeriodicity];	// 6 terms can be stored
	double		_Vs_kJ[MaxPeriodicity];		// 6 V terms
        double          _PhaseRads[MaxPeriodicity]; 	// 6 phase terms

	void	setTypes( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4);
	int	maxPeriodicity() const { return MaxPeriodicity; };
        bool	hasPeriodicity(int period) const;
        double  getV_kJ(int period) const;
        void    setV_kJ(int period, double v);
	double	getV_kCal(int period) const;
	void	setV_kCal(int period, double v);

        double  getPhaseRad(int period) const;
        void    setPhaseRad(int period, double pr);

        void    mergeWith(FFPtor_sp ptor);

	virtual	ParameterType	type() { return ptor;};
	virtual	string		levelDescription();
	DEFAULT_CTOR_DTOR(FFPtor_O);
    };




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
    SMART(FFPtorDb);
    class FFPtorDb_O : public FFBaseDb_O
    {
	LISP_BASE1(FFBaseDb_O);
	LISP_CLASS(chem,ChemPkg,FFPtorDb_O,"FFPtorDb");

    public:
#ifdef XML_ARCHIVE
	void	archiveBase(core::ArchiveP node);
#endif

    public:
	typedef	gctools::Vec0<FFPtor_sp>::iterator	FFPtor_spIterator;
        gctools::Vec0<FFPtor_sp>        _Terms;
	adapt::SymbolMap<FFPtor_O>    _Lookup;

	FFPtor_spIterator	begin() { return this->_Terms.begin();};
	FFPtor_spIterator	end() { return this->_Terms.end();};

	void	setTermFormat( core::VectorStrings vterms );

	void	add( FFPtor_sp ptor );

	FFPtor_sp findExactTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );
	bool    hasExactTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );

	//! Look for exact term then general one
        FFPtor_sp findBestTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );
        bool    hasBestTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );

        void    cantFind(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 );

	//! Dump all ptors to stdout that match these types
//	void	list( string a1, string a2, string a3, string a4 );


	void		initialize();


	DEFAULT_CTOR_DTOR(FFPtorDb_O);
    };


};

TRANSLATE(chem::FFPtor_O);
TRANSLATE(chem::FFPtorDb_O);
#endif
