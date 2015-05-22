       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	macroModel.h
 *
 *	Functions that read and write macroModel files.
 *
 */

#ifndef	MACROMODEL_H
#define	MACROMODEL_H

#include <iostream>
#include <string>
#include "core/common.h"
#include "matter.h"
#include "molecule.h"
#include "residue.h"
#include "atom.h"


#include "chemPackage.h"


namespace chem {


    FORWARD(MacroModelFile);
    class MacroModelFile_O : public core::T_O
    {

	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,MacroModelFile_O,"MacroModelFile");
public:

	void	openRead(const string& fn );
	void	openWrite(const string& fn );
	void	close();

	bool		atEnd();

	Molecule_sp	readMolecule();
	void		writeAggregate( Aggregate_sp );
private:
	int		lines;
	FILE*		file;

	DEFAULT_CTOR_DTOR(MacroModelFile_O);
};

};




#endif
