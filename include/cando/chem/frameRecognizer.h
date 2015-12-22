       
       

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
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(chem,ChemPkg,FrameRecognizer_O,"FrameRecognizer");
#if INIT_TO_FACTORIES
 public:
    static FrameRecognizer_sp make(core::Symbol_sp name, const string& smarts, core::Symbol_sp groupName);
#else
    DECLARE_INIT();
#endif
public:
	void	initialize();
//	void	archive(core::ArchiveP node);

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
	ChemInfo_sp		_ChemInfo;
public:

public: // Lisp creation
	static core::T_sp prim_LispCreate(core::Function_sp e, core::List_sp args);


public:	// ----------------------------------------------------------


    void	setGroupName(core::Symbol_sp fn);
    core::Symbol_sp	getGroupName();

    
    void	setRecognizerName(core::Symbol_sp fn);
    core::Symbol_sp	getRecognizerName();

CL_NAME("getSmarts");
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
