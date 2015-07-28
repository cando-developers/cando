       
       
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
#include <clasp/core/common.h>
#include <cando/chem/matter.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>


#include <cando/chem/chemPackage.h>


namespace chem {


    FORWARD(MacroModelFile);
    class MacroModelFile_O : public core::CxxObject_O
    {

	LISP_BASE1(core::CxxObject_O);
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
