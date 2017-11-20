/*
    File: candoDatabaseDependent.h
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
       
       
#ifndef	CandoDatabaseDependent_H //[
#define CandoDatabaseDependent_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>



#include <cando/chem/chemPackage.h>


namespace chem {

SMART(CandoDatabase);
SMART(Constitution);

#if 0
SMART(CandoDatabaseDependent);
__ BEGIN_CLASS_DEFINITION(CandoDatabase_ODependent,CandoDatabaseDependent,O EnvironmentDependent) // {
public:
	void initialize();
public:
//	void	archiveBase(core::ArchiveP node);
//	bool	loadFinalize(core::ArchiveP node);
private:
//	CandoDatabase_sp	_WeakCandoDatabase;
public:


//    virtual void	setCandoDatabase(CandoDatabase_sp bdb);

		/*! Return true if this name can not be expanded.
		 * For instance monomer groups can be expanded to every monomer within them.
		 */
    virtual bool	isTerminalName() { return false; };
    		/*! Get the name of this object
		 */
    virtual string	getName()	{_OF();SUBCLASS_MUST_IMPLEMENT();};

    		/*! Return a list of names that is expanded recursively as much as possible
		 */
    virtual adapt::SymbolSet_sp	expandedNameSet() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
    		/*! Return a list of names that is expanded recursively but terminate with
		 * representative names.
		 */
    virtual RepresentativeList_sp expandedRepresentativeList() const {_OF(); SUBCLASS_MUST_IMPLEMENT(); };

		/*! If this object is a Stereoisomer then it has a Constitution otherwise it doesn't */
    virtual bool hasConstitution() { return false; };
    		/*! Return the constitution */
    virtual Constitution_sp constitution();


	CandoDatabase_sp getCandoDatabase();
	bool hasCandoDatabase();


	CandoDatabase_ODependent( const CandoDatabase_ODependent& ss ); //!< Copy constructor


__END_CLASS_DEFINITION(CandoDatabase_ODependent) //}
#endif

};
#endif //]
