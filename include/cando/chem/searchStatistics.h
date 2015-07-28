#ifndef	SearchStatistics_H //[
#define SearchStatistics_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/numerics.h>


#include <cando/chem/chemPackage.h>


namespace chem {



    SMART(SearchStatistics );
    class SearchStatistics_O : public core::CxxObject_O
    {
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(chem,ChemPkg,SearchStatistics_O,"SearchStatistics");
#if INIT_TO_FACTORIES
    public:
	static SearchStatistics_sp make();
#else
	DECLARE_INIT();
#endif
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;
	string description() const;

private: // instance variables
	Bignum		_Builds;
	Bignum		_Rejects;
	Bignum		_Hits;
	uint			_ReportFrequency;

public:	// Creation class functions

public:
	void setReportFrequency(uint num);

	void reset();

	void incrementBuilds();
	Bignum getBuilds();

	void incrementRejects();
	Bignum getRejects();

	void incrementHits();
	Bignum getHits();

	SearchStatistics_O( const SearchStatistics_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(SearchStatistics_O);
    };


};
TRANSLATE(chem::SearchStatistics_O);
#endif //]
