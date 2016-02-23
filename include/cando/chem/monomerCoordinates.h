/*
    File: monomerCoordinates.h
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
       
       



#ifndef	MonomerCoordinates_H //[
#define MonomerCoordinates_H


#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/conformationExplorer.h>



#include <cando/chem/chemPackage.h>


namespace chem {

    FORWARD(Rotamer);
SMART(Topology);

SMART(MonomerCoordinates);
class MonomerCoordinates_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,MonomerCoordinates_O,"MonomerCoordinates",core::CxxObject_O);
#if INIT_TO_FACTORIES
 public:
    static MonomerCoordinates_sp make(MonomerContext_sp context);
#else
    DECLARE_INIT();
#endif
public:
	void	initialize();
//	void	archiveBase(core::ArchiveP node);

public:
    typedef	gctools::Vec0<Rotamer_sp>	Rotamers;
private:
	string				_Comment;
	MonomerContext_sp		_Context;
	string				_TrainerName;
	Rotamers			_Rotamers;
public:
	void _defineFromConformationExplorerOrDebug( Topology_sp topology,
				ConformationExplorer_sp structureList,
				Monomer_sp  focusMonomer,
				bool	   debug,
				uint maximumNumberOfConformations );
public:

CL_NAME("setComment");
CL_DEFMETHOD 	void setComment(const string& context) { this->_Comment = context; };
CL_NAME("getComment");
CL_DEFMETHOD 	string getComment() { return this->_Comment; };

	void setContext(MonomerContext_sp context) { this->_Context = context; };
CL_NAME("getContext");
CL_DEFMETHOD 	MonomerContext_sp getContext() { return this->_Context; };

	void defineFromConformationExplorer( Topology_sp topology,
				ConformationExplorer_sp structureList,
				Monomer_sp  focusMonomer,
				uint maximumNumberOfConformations );

	/*! Test if the extraction process will work with the current Aggregate
	 * if it fails then there is a problem in the MetaDatabase that needs to
	 * be fixed before we go forward
	 * Will return a Cons with two KeyedObjects (sawProblem: (bool))
	 * and (problems: (StringStream))
	 */

	 core::List_sp testExtraction( Topology_sp topology, 
			Aggregate_sp aggregate, 
			Monomer_sp focusMonomer );
	
	DEFAULT_CTOR_DTOR(MonomerCoordinates_O);
};




};
TRANSLATE(chem::MonomerCoordinates_O);
#endif //]
