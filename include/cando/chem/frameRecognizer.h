/*
    File: frameRecognizer.h
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
       
       

#ifndef	FrameRecognizer_H //[
#define	FrameRecognizer_H


#include <clasp/core/common.h>
#include <string>
////#include "core/archiveNode.h"
#include <cando/chem/chemInfo.h>

#include <cando/chem/chemPackage.h>


namespace chem {

class	BoundFrame_O;
    typedef	gctools::smart_ptr<BoundFrame_O>	BoundFrame_sp;

SMART(FrameRecognizer);
class FrameRecognizer_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,FrameRecognizer_O,"FrameRecognizer",core::CxxObject_O);
#if INIT_TO_FACTORIES
 public:
    static FrameRecognizer_sp make(core::Symbol_sp name, const string& smarts, core::Symbol_sp groupName);
#else
    DECLARE_INIT();
#endif
public:
	void	initialize();

private: // ---------------------------------------------------------
	core::Symbol_sp			_Name;
		/*! The group name is used to group together FrameRecognizers
		 * that are similar to each other
		 * eg: a frame recognizer that recognizes a primary amide and a secondary
		 * amide.
		 * If the GroupName is empty then the _Name is returned.
		 */
	core::Symbol_sp		_GroupName;
	string			_Smarts;

private: // Do not archive
  Root_sp		_Root;
public:

public: // Lisp creation
	static core::T_sp prim_LispCreate(core::Function_sp e, core::List_sp args);


public:	// ----------------------------------------------------------


    void	setGroupName(core::Symbol_sp fn);
    core::Symbol_sp	getGroupName();

    
    void	setRecognizerName(core::Symbol_sp fn);
    core::Symbol_sp	getRecognizerName();

CL_LISPIFY_NAME("getSmarts");
CL_DEFMETHOD     string	getSmarts()	{ return this->_Smarts;};
    void	compileSmarts(const string& fn);

	/*! Given the O atom, return true if this FrameRecognizer recognizes
	 * the neighboring atoms and set the temporaryPAtomMatch and temporaryQAtomMatch
	 * variables
	 */
    bool recognizes( Atom_sp atomO );

    string description() const;
	/*! Return the ChemInfoMatch object for the match
	 */
    ChemInfoMatch_sp getMatch();
    	/*! Return the depth of the longest pattern match
	 */
    uint depth();


    DEFAULT_CTOR_DTOR(FrameRecognizer_O);
};



};


TRANSLATE(chem::FrameRecognizer_O);
#endif
