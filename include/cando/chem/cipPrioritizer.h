#ifndef	CipPrioritizer_H //[
#define CipPrioritizer_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>


namespace chem {

SMART(Matter);
SMART(Atom);

SMART(CipPrioritizer );
class CipPrioritizer_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,CipPrioritizer_O,"CipPrioritizer",core::CxxObject_O);
private:
	// instance variables
    std::vector<int>		_p;
    std::vector< std::vector<int> > 	_s;
public:
	void	initialize();
public:
    static void assignPriorities(Matter_sp molOrAgg);
public:
	int	getP(Atom_sp a);
    std::vector<int>& getS(Atom_sp a);
     void assignCahnIngoldPrelogPriorityToAtomsRelativePriority(Matter_sp molOrAgg);
public:

	void defineStereochemicalConfigurationsForAllAtoms(Matter_sp molOrAgg);
	CipPrioritizer_O( const CipPrioritizer_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(CipPrioritizer_O);
};


};
TRANSLATE(chem::CipPrioritizer_O);
#endif //]
