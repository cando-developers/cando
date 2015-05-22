       
       
       
/*!
__BEGIN_DOC( msmarts, MSMARTS chemical pattern matching)

Based on SMARTS documentation at 
\begin{verbatim}
http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
\end{verbatim}

MSMARTS is similar to SMARTS with the following differences.
\begin{itemize}
	\item MSMARTS supports atom tags: numerical labels that can be attached to atoms as an MSMARTS substructure is matched to a molecule. The tagged atoms can then be referenced after the substructure is matched.

	For example: the MSMARTS string "[N\&H1]1C2(=O3)" will recognize a secondary amide and the amide nitrogen, carbonyl carbon and carbonyl oxygen can be obtained using the tags "1", "2" and "3" after a successful match.
	\item The syntax for identifying rings is different. Rings are recognized with strings like: "C1CCC[C\&?1]".  The first "1" assigns a tag "1" to the first carbon, The "[C\&?1]" atom tests if the atom is carbon and has the tag "1".
\end{itemize}

   Substructure searching, the process of finding a particular pattern (subgraph) in a molecule (graph), is one of the most important tasks for computers in chemistry. It is used in virtually every application that employs a digital representation of a molecule, including depiction (to highlight a particular functional group), drug design (searching a database for similar structures and activity), analytical chemistry (looking for previously-characterized structures and comparing their data to that of an unknown), and a host of other problems.

   MSMARTS expressions allow a chemist to specify substructures using rules that 
   are straightforward extensions of SMILES. For example: to search a database 
   for phenol-containing structures, one would use the SMARTS string 
   [OH]c1cccc[c\&?1], which is similar to SMILES (Note: the [c\&?1] atom 
   primative test is used to identify rings in MSMARTS). In fact, almost all 
   SMILES specifications are valid SMARTS targets. Using SMARTS, flexible 
   and efficient substructure-search specifications can be made in terms 
   that are meaningful to chemists.

   In the SMILES language, there are two fundamental types of symbols: atoms and bonds. Using these SMILES symbols, once can specify a molecule's graph (its "nodes" and "edges") and assign "labels" to the components of the graph (that is, say what type of atom each node represents, and what type of bond each edge represents).

   The same is true in SMARTS: One uses atomic and bond symbols to specify a graph. However, in SMARTS the labels for the graph's nodes and edges (its "atoms" and "bonds") are extended to include "logical operators" and special atomic and bond symbols; these allow SMARTS atoms and bonds to be more general. For example, the SMARTS atomic symbol [C,N] is an atom that can be aliphatic C or aliphatic N; the SMARTS bond symbol ~ (tilde) matches any bond.

Below is example code that uses SMARTS to find every amide bond in a molecule:
\begin{verbatim}
#
# Define a ChemInfo object that can carry out
# substructure searches
#
( amideSmarts = [ new ChemInfo] )
#
# Compile a substructure pattern using SMARTS code
#
( amideSmarts compileSmarts "N1~C2=O3" )
# 
# Load a molecule
#
( p53 = [ loadMol2 "p53.mol2" ] )
#
# Iterate through every atom and if it matches
# the substructure search then extract the tagged
# atoms and print their names
#
(for a (atomsAsCons p53)
    (if (amideSmarts matches a)
        (progn
		(print (% "-----Matching atom: %s" (getName a)))
		(let ((coAtom (getAtomWithTag amideSmarts "2"))
		      (oAtom (getAtomWithTag amideSmarts "3"))) 
		   (print (% "    Carbonyl carbon: %s" (getName coAtom)))
		   (print (% "    Carbonyl oxygen: %s" (getName oAtom)))))))
\end{verbatim}
__END_DOC
*/



%{
#define DEBUG_LEVEL_FULL

#include<iostream>
#include<string>
#include<vector>
#include <istream>

#include "core//common.h"
#include "aggregate.h"
#include "molecule.h"
#include "residue.h"
#include "atom.h"

#include "chemInfo.h"
#include "hold.h"


using namespace std;
using namespace chem;

// This defines the stream that we are reading from

#include "msmarts_TypeParser.h"
#include "msmarts_ParserParam.h"
#include "msmarts_LocationType.h"

#define	LEXDEBUG	1

#ifdef	LEXDEBUG
#define	LEXPRINT(x) {lisp_LOG(BF( "Token: %s") %(x));}
#define	LEXDPRINT(x) {lisp_LOG(x)};
#else
#define	LEXPRINT(x)
#define	LEXDPRINT(x)
#endif
#define	MAX_CVAL	16

#define MP() (data)

void msmarts_error(YYLTYPE* yyloc, msmarts_SParserParam* data, const char* message );
int	msmarts_lex(YYSTYPE* yylval, YYLTYPE* yylloc, msmarts_SParserParam* data);

%}

%name-prefix="msmarts_"
%define api.pure
%locations
%lex-param {msmarts_SParserParam* data}
%parse-param {msmarts_SParserParam* data}


%token APEndOfLine
%token <benum> APBond
%token APAtomicNumber
%token APChiralityAntiClockwise
%token APChiralityClockwise
%token APCloseBracket
%token APCloseParenthesis
%token APConnectivity
%token APDegree
%token APRingTest
%token APResidueTest
//%token APDirectionalSingleDown
//%token APDirectionalSingleUp
%token APDollarSign
%token APDoubleBond
%token APAliphatic
%token APAromatic
%token <eval> APOrganicElement
%token <eval> APAromaticOrganicElement
%token <eval> APInorganicElement
%token APElectronegativeElement
%token APError
%token APGroupNumber  
%token APHeavyAtomTotalBond
%token APImplicitHCount
%token AM1_BCC_x
%token AM1_BCC_y
%token APLonePair
%token APNegativeCharge
%token APNegativeFormalCharge   
%token <cval> APNumber
%token <cval> APLambda
%token APOpenBracket
%token APOpenParenthesis
%token APOperatorAndHigh
%token APOperatorAndLow
%token APOperatorNot
%token APOperatorOr
%token APPiBondOrbital
%token APAromaticPiElectron
%token APPositiveCharge
%token APPositiveFormalCharge
%token APRingMembershipCount
%token APRingSize
//%token APSingleBond
//%token APSingleDownOrUnspecified
//%token APSingleUpOrUnspecified
%token APTotalHCount
%token APTotalBondNumber		//Explicit + Implicit bonds 
%token APTotalExplicitBondNumber	//Explicit bonds only
%token APTransitionMetal
//%token APTripleBond
%token APValence
%token APWildCard

%left	APOperatorAndLow
%left	APOperatorOr
%left	APOperatorAndHigh
%left	APOperatorNot



//%type <smartsRoot>		input
%type <smartsRoot>		chemInfo
%type <root>			recursiveChemInfo
%type <ival>		intNumber
%type <bondListMatchNode>	chain branch
%type <atomOrBondMatchNode>	atomTest atomPrimativeTest
%type <bondTest>		bondAtomTest
%type <atomOrBondMatchNode>	acyclicAtomTest
%type <logical>			logOp

//%destructor	{ if ($$!=NULL) delete ($$);} 	input
%destructor	{ if ($$!=NULL) delete ($$);} 	chemInfo
%destructor { if ($$!=NULL) delete ($$); } 	recursiveChemInfo
%destructor { if ($$!=NULL) delete ($$); } 	chain
%destructor { if ($$!=NULL) delete ($$); } 	branch
%destructor { if ($$!=NULL) delete ($$); } 	atomTest
%destructor { if ($$!=NULL) delete ($$); } 	atomPrimativeTest
%destructor { if ($$!=NULL) delete ($$); } 	bondAtomTest
%destructor { if ($$!=NULL) delete ($$); } 	acyclicAtomTest
%destructor { if ($$!=NULL) delete ($$); } 	logOp
%%


input:  chemInfo	
	{ _lisp_BLOCK_TRACE("input: chemInfo");
	    data->expression = $1->_obj; 
	};

chemInfo:	/* empty */ { $$ = new Hold<chem::SmartsRoot_O>(); }
	    | atomTest 	
	    	{ _lisp_BLOCK_TRACE("chemInfo: atomTest");
		    $$ = new Hold<chem::SmartsRoot_O>(chem::SmartsRoot_O::create($1->_obj)); 
		};
	    | atomTest chain 
	    	{ _lisp_BLOCK_TRACE("chemInfo:atomTest chain");
		    $$ = new Hold<chem::SmartsRoot_O>(chem::SmartsRoot_O::create($1->_obj,$2->_obj)); 
		}
	    ;


chain: bondAtomTest chain 
    { _lisp_BLOCK_TRACE("chain:bondAtomTest chain");
	$$ = new Hold<chem::BondListMatchNode_O>(chem::Chain_O::create( $1->_obj, $2->_obj)); 
    }
| branch chain	   
    { _lisp_BLOCK_TRACE("chain: branch chain");
	$$ = new Hold<chem::BondListMatchNode_O>(chem::Branch_O::create( $1->_obj, $2->_obj)); 
    }
| bondAtomTest	   
    { _lisp_BLOCK_TRACE("chain: bondAtomTest");
	$$ = new Hold<chem::BondListMatchNode_O>(chem::Chain_O::create( $1->_obj)); 
    }
| branch  
    { _lisp_BLOCK_TRACE("chain: branch");
	$$ = new Hold<chem::BondListMatchNode_O>(chem::Branch_O::create( $1->_obj)); 
    } /* handle C(C=O) */
;

bondAtomTest:	APBond atomTest 
			{ _lisp_BLOCK_TRACE("bondAtomTest:APBond atomTest");
			    $$ = new Hold<chem::BondTest_O>(chem::BondTest_O::create( $1, $2->_obj)); 
			}
	    |	atomTest 
	    		{ _lisp_BLOCK_TRACE("bondAtomTest: atomTest");
			    $$ = new Hold<chem::BondTest_O>(chem::BondTest_O::create( chem::SABSingleOrAromaticBond, $1->_obj)); 
			} 
	    ;
    //	|	APSingleBond atomTest { $$ = new Holder(chem::BondTest_O::create( chem::SABSingleBond, $2->_obj)); }
    //	|	APDoubleBond atomTest { $$ = new Holder(chem::BondTest_O::create( chem::SABDoubleBond, $2)->_obj); }
    // 	|	APTripleBond atomTest { $$ = new Holder(chem::BondTest_O::create( chem::SABTripleBond, $2)->_obj); }
    //	|	APDirectionalSingleUp atomTest { 
    //			$$ = chem::BondTest_O::create( chem::SABDirectionalSingleUp, $2);  }
    //	|	APSingleUpOrUnspecified atomTest { 
    //			$$ = chem::BondTest_O::create(chem::SABDirectionalSingleUpOrUnspecified, $2);}
    //	|       APDirectionalSingleDown atomTest { 
    //			$$ = chem::BondTest_O::create( chem::SABDirectionalSingleDown, $2); } 
    //	|	APSingleDownOrUnspecified atomTest { 
    //			$$ = chem::BondTest_O::create( chem::SABDirectionalSingleDownOrUnspecified,$2);}
    //	|	APAromaticBond atomTest { 
    //			$$ = chem::BondTest_O::create( chem::SABAromaticBond, $2);}



acyclicAtomTest: APOpenBracket logOp APCloseBracket 
	{ _lisp_BLOCK_TRACE("acyclicAtomTest: APOpenBracket logOp APCloseBracket");
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>($2->_obj);
	} 
| APOrganicElement 
 	{ _lisp_BLOCK_TRACE("acyclicAtomTest: APOrganicElement");
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>((chem::AtomTest_O::create(chem::SAPElement, _lisp->internKeyword($1)))); 
	} 

| APAromaticOrganicElement 
 	{ _lisp_BLOCK_TRACE("acyclicAtomTest: APAromaticOrganicElement");
	    string elementName = $1;
	    elementName[0] = toupper(elementName[0]);
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>((chem::AtomTest_O::create(chem::SAPAromaticElement, _lisp->internKeyword(elementName)))); 
	} 

| APAliphatic
	{
	    $$=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPAliphatic));
	}

| APAromatic
	{
	    $$=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPAromatic));
	}


| recursiveChemInfo 
	    		{ _lisp_BLOCK_TRACE("acyclicAtomTest: recursiveChemInfo");
	    		    $$ = new Hold<chem::AtomOrBondMatchNode_O>($1->_obj); 
			}
	    ;

atomTest: acyclicAtomTest APNumber 
    { _lisp_BLOCK_TRACEF(BF("atomTest: acyclicAtomTest APNumber   // atomTag==(%s)  ") % ($2));
	$$=new Hold<chem::AtomOrBondMatchNode_O>(chem::TagSet_O::create($1->_obj,$2));
    }
| acyclicAtomTest 
    { _lisp_BLOCK_TRACE("atomTest: acyclicAtomTest");
	$$ = new Hold<chem::AtomOrBondMatchNode_O>($1->_obj); 
    };

branch:	APOpenParenthesis chain APCloseParenthesis { $$ = new Hold<chem::BondListMatchNode_O>($2->_obj); };


 /* 
__BEGIN_DOC( msmarts.atomics, Atomic Primitives)
   SMARTS provides a number of primitive symbols describing atomic properties
   beyond those used in SMILES (atomic symbol, charge, and isotopic
   specifications). The following tables list the atomic primitives
   used in SMARTS (all SMILES atomic symbols are also legal).
   In these tables <n> stands for a digit, <c> for chiral class.

   Note that atomic primitive H can have two meanings,
   implying a property or the element itself. [H] means
   hydrogen atom. [*H2] means any atom with exactly two hydrogens attached
   
   	\begin{tabular}{| l | l | l | l |}
	\hline
	Symbol & Symbol name & Atomic property requirements & Default \\ \hline
__END_DOC
 */

 /*
__APPEND_DOC(msmarts.atomics)
	* & wildcard & any atom & (no default) \\ \hline
__END_DOC
 */
atomPrimativeTest :	APWildCard 
	{ 
		$$ = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPWildCard)); 
	}

/*
__APPEND_DOC(msmarts.atomics)
	D\emph{n} & APDegree & explicit connections & exactly one \\ \hline
__END_DOC
*/
| APDegree 
    	{ 
		$$ = new Hold<chem::AtomOrBondMatchNode_O>( 
				chem::AtomTest_O::create( chem::SAPDegree, 1)); 
	}

| APDegree intNumber 
	{ 
		$$ = new Hold<chem::AtomOrBondMatchNode_O>( 
			chem::AtomTest_O::create( chem::SAPDegree,$2)); 
	}


/*
__APPEND_DOC(msmarts.atomics)
	x & AM1_BCC_x & AM1-BCC-x aromatic atom test: J Comput Chem 23: 1623–1641, 2002 \\ \hline
__END_DOC
*/
| AM1_BCC_x
	{ 
		$$ = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create( chem::SAPAM1_BCC_x, 0));
	}

/*
__APPEND_DOC(msmarts.atomics)
	x & AM1_BCC_y & AM1-BCC-y aromatic atom test: J Comput Chem 23: 1623–1641, 2002 \\ \hline
__END_DOC
*/
| AM1_BCC_y
	{ 
		$$ = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create( chem::SAPAM1_BCC_y, 0));
	}


/*
__APPEND_DOC(msmarts.atomics)
	H\emph{n} & APTotalHCount& \emph{n} attached hydrogens & exactly one \\ \hline
__END_DOC
*/
| APTotalHCount 
	{ 
		$$ = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create( chem::SAPTotalHCount, 1));
	}
| APTotalHCount intNumber 
	{ 
		$$ = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create( chem::SAPTotalHCount, $2));
	}

/*
__APPEND_DOC(msmarts.atomics)
	h\emph{n} & APImplicitHCount& \emph{n} implicit attached hydrogens & at least one\\ \hline
__END_DOC
*/
| APImplicitHCount 
	{
		$$ = new Hold<chem::AtomOrBondMatchNode_O>( 
			chem::AtomTest_O::create( chem::SAPImplicitHCount, 1)); 
	}
| APImplicitHCount intNumber 
	{ 
		$$ = new Hold<chem::AtomOrBondMatchNode_O>(
					chem::AtomTest_O::create( 
						chem::SAPImplicitHCount, $2));
	}

/*
__APPEND_DOC(msmarts.atomics)
	?\emph{n} & APRingTest& Atom is matched to atom tagged with \emph{n}& (no default)\\ \hline
__END_DOC
*/
| APRingTest APNumber
	{
		$$=new Hold<chem::AtomOrBondMatchNode_O>(
					chem::AtomTest_O::create(chem::SAPRingTest,$2));
	}


/*
__APPEND_DOC(msmarts.atoms)
	A & APAliphatic & True if not aromatic
__END_DOC
*/
| APAliphatic
	{
	    $$=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPAliphatic));
	}

/*
__APPEND_DOC(msmarts.atoms)
	a & APAromatic & True if aromatic
__END_DOC
*/
| APAromatic
	{
	    $$=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPAromatic));
	}

/*
__APPEND_DOC(msmarts.atoms)
	<xxx> & APLambda & Call lambda with name xxx with atom as arg \hline
__END_DOC
*/

| APLambda
	{
	    string symStr($1);
	    core::Symbol_sp sym = _lisp->intern(symStr);
	    $$=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPLambda,sym));
	}


/*
__APPEND_DOC(msmarts.atomics)
	U\emph{n} & APResidueTest& Atom must be in same residue as atom tagged \emph{n}& (no default)\\ \hline
__END_DOC
*/
| APResidueTest APNumber
	{
		$$=new Hold<chem::AtomOrBondMatchNode_O>(
					chem::AtomTest_O::create(chem::SAPResidueTest,$2));
	}


/*
__APPEND_DOC(msmarts.atomics)
	R\emph{n} & APRingMemberCount & is in \emph{n} SSSR rings \par(WORKS?)& any ring atom\\ \hline
__END_DOC
*/
| APRingMembershipCount 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPRingMembershipCount, 1)); 
	}
| APRingMembershipCount intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPRingMembershipCount, $2));
	}


/*
__APPEND_DOC(msmarts.atomics)
	r\emph{n} & APRingSize & is in smallest SSSR size \emph{n} & any ring atom\\ \hline
__END_DOC
*/
| APRingSize 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPRingSize, 1)); 
	}	
| APRingSize intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPRingSize, $2));
	}
/*
__APPEND_DOC(msmarts.atomics)
	v\emph{n} & APValence & total bond order \emph{n} & exactly 1 \\ \hline
__END_DOC
*/
| APValence 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPValence, 1)); 
	}
| APValence intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPValence, $2));
	}

/*
__APPEND_DOC(msmarts.atomics)
	X\emph{n} & APConnectivity& \emph{n} total connections & exactly 1 \\ \hline
__END_DOC
*/
| APConnectivity 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPConnectivity, 1)); 
	}
| APConnectivity intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPConnectivity, $2));
	}


| APGroupNumber intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPGroupNumber, $2)); 
	} 
| APElectronegativeElement 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPElectronegativeElement)); 
	} 
| APHeavyAtomTotalBond 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPHeavyAtomTotalBond, 1)); 
	}
| APHeavyAtomTotalBond intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPHeavyAtomTotalBond, $2)); 
	} 
/*| APNegativeCharge APWildcardOrAtomicMass intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeFormalCharge, $3 )); 
	}  */
/*
__APPEND_DOC(msmarts.atomics)
	-\emph{n} & APNegativeCharge & -\emph{n} charge & exactly -1 \\ \hline
__END_DOC
*/
| APNegativeCharge intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeCharge, $2));
	}
| APNegativeCharge 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeCharge, -1));
	}

/*
__APPEND_DOC(msmarts.atomics)
	-- & APNegativeCharge 2x & -2 charge & exactly -2 \\ \hline
__END_DOC
*/
| APNegativeCharge APNegativeCharge
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeCharge, -2));
	}
/*
__APPEND_DOC(msmarts.atomics)
	--- & APNegativeCharge 3x & -3 charge & exactly -3 \\ \hline
__END_DOC
*/
| APNegativeCharge APNegativeCharge APNegativeCharge 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeCharge, -3));
	}
/*| APPositiveCharge APWildcardOrAtomicMass intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveFormalCharge, $3 )); 
	}  */



/*
__APPEND_DOC(msmarts.atomics)
	+\emph{n} & APPositiveCharge & +\emph{n} charge & exactly +1 \\ \hline
__END_DOC
*/
| APPositiveCharge intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveCharge, $2));
	}
| APPositiveCharge 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveCharge, 1));
	}
/*
__APPEND_DOC(msmarts.atomics)
	++ & APPositiveCharge 2x & +2 charge & exactly +2 \\ \hline
__END_DOC
*/
| APPositiveCharge APPositiveCharge
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveCharge, 2));
	}
/*
__APPEND_DOC(msmarts.atomics)
	+++ & APPositiveCharge 3x & +3 charge & exactly +3 \\ \hline
__END_DOC
*/
| APPositiveCharge APPositiveCharge APPositiveCharge 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveCharge, 3));
	}
| APTransitionMetal intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPTransitionMetal, $2)); 
	}
| APAromaticPiElectron intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPAromaticPiElectron , $2));
	}
| APPiBondOrbital 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPiBondOrbital)); 
	}
/*
__APPEND_DOC(msmarts.atomics)
	\#\emph{n} & APAtomicNumber & atomic number \emph{n} & (no default) \\ \hline
__END_DOC
*/
| APAtomicNumber intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPAtomicNumber, $2));
	}
| APChiralityAntiClockwise 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPChiralityAntiClockwise));
	}
| APChiralityClockwise 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPChiralityClockwise));
	}
/*
__APPEND_DOC(msmarts.atomics)
	\emph{n} & APAtomicMass & atomic mass \emph{n} & (no default) \\ \hline
__END_DOC
*/
| intNumber 
	{
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPAtomicMass, $1));
	}
| APOrganicElement 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPElement, _lisp->internKeyword($1) )); 
	}
| APAromaticOrganicElement 
	{ _lisp_BLOCK_TRACE("atomPrimitiveTest: APAromaticOrganicElement");
	    string elementName = $1;
	    elementName[0] = toupper(elementName[0]);
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPAromaticElement, _lisp->internKeyword(elementName) )); 
	}
| APInorganicElement 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPElement, _lisp->internKeyword($1) )); 
	}
| APLonePair 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPLonePair, 1)); 
	}
| APLonePair intNumber 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPLonePair, $2));
	}
/*
__APPEND_DOC(msmarts.atomics)
	\$(\emph{\tiny MSMARTS}) & recursive MSMARTS & match recursive MSMARTS & (no default) \\ \hline
__END_DOC
*/
| recursiveChemInfo 
	{ 
	    $$ = new Hold<chem::AtomOrBondMatchNode_O>($1->_obj); 
	}
    ;
/*
__APPEND_DOC(msmarts.atomics)
\end{tabular}

Some of these have not been debugged. Test before you trust them.

Examples: \par
\begin{tabular}{| l | l |}
\hline
  [CH2] & aliphatic carbon with two hydrogens (methylene carbon) \\ \hline
  [!C;R] & ( NOT aliphatic carbon ) AND in ring \\ \hline
  [!C;!R0] & same as above ("!R0" means not in zero rings) \\ \hline
  [n;H1] & H-pyrrole nitrogen \\ \hline
  [n\&H1] & same as above \\ \hline
  [nH1] & same as above \\ \hline
  [c,n\&H1] & any arom carbon OR H-pyrrole nitrogen \\ \hline
  [X3\&H0] & atom with 3 total bonds and no H's \\ \hline
  [c,n;H1] & (arom carbon OR arom nitrogen)  and exactly one H \\ \hline
  [Cl] & any chlorine atom \\ \hline
  [35*] & any atom of mass 35 \\ \hline
  [35Cl] & chlorine atom of mass 35 \\ \hline
  [F,Cl,Br,I] & the 1st four halogens. \\ \hline
\end{tabular}
__END_DOC
*/

/*
__BEGIN_DOC(msmarts.recursive,Recursive MSMARTS)

     Any MSMARTS expression may be used to define an atomic environment by writing a SMARTS starting with the atom of interest in this form:
     \$(\emph{MSMARTS})
     Such definitions may be considered atomic properties. These expressions can be used in same manner as other atomic primitives (also, they can be nested). Recursive SMARTS expressions are used in the following manner:

	\begin{tabular}{| l | l |}
	\hline
	*C & atom connected to methyl or methylene carbon \\ \hline
	*CC & atom connected to ethyl carbon \\ \hline
	[\$(*C);\$(*CC)] & Atom in both above environments (matches CCC) \\ \hline
	\end{tabular}

The additional power of such expressions is illustrated by the following
example which derives an expression for methyl carbons which are ortho to
oxygen and meta to a nitrogen on an aromatic ring.


\begin{tabular}{| l | l | }
\hline
CaaO & C ortho to O \\ \hline
CaaaN & C meta to N \\ \hline
Caa(O)aN & C ortho to O and meta to N (but 2O,3N only) \\ \hline
Ca(aO)aaN & C ortho to O and meta to N (but 2O,5N only) \\ \hline
C[\$( aaO);\$( aaaN)] & C ortho to O and meta to N (all cases) \\ \hline
\end{tabular}
__END_DOC
 */


 /*!
__BEGIN_DOC(msmarts.logical,Logical Operators)

 Atom and bond primitive specifications may be combined to form expressions by using logical operators. In the following table, e is an atom or bond SMARTS expression (which may be a primitive). The logical operators are listed in order of decreasing precedence (high precedence operators are evaluated first).

All atomic expressions which are not simple primitives must be enclosed in brackets. The default operation is \& (high precedence "and"), i.e., two adjacent primitives without an intervening logical operator must both be true for the expression (or subexpression) to be true.

The ability to form expressions gives the SMARTS user a great deal of power to specify exactly what is desired. The two forms of the AND operator are used in SMARTS instead of grouping operators.

\begin{tabular}{|l|l|l|}
\hline
\textbf{Symbol} & \textbf{Expression} & \textbf{Meaning} \\ \hline
__END_DOC
*/

logOp:	atomPrimativeTest 	
	{ 
	    $$ = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logIdentity, $1->_obj )); 
	}
/*!
__APPEND_DOC(msmarts.logical)
exclamation & !e1 & not e1 \\ \hline
__END_DOC
*/
| APOperatorNot logOp 
	{ 
	    $$ = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logNot, $2->_obj )); 
	}
/*
__APPEND_DOC(msmarts.logical)
ampersand & e1\&e2 & e1 and e2 (high precedence) \\ \hline
__END_DOC
*/
| logOp APOperatorAndHigh logOp 
	{ 
	    $$ = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logHighPrecedenceAnd, $1->_obj, $3->_obj )); 
	}
/*
__APPEND_DOC(msmarts.logical)
comma & e1,e2  & e1 or e2 \\ \hline
__END_DOC
*/
| logOp APOperatorOr logOp 
	{ 
	    $$ = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logOr, $1->_obj, $3->_obj )); 
	}
/*
__APPEND_DOC(msmarts.logical)
semicolon & e1;e2  & e1 and e2 (low precedence) \\ \hline
\end{tabular}
__END_DOC
*/
| logOp APOperatorAndLow logOp 
	{ 
	    $$ = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logLowPrecedenceAnd, $1->_obj, $3->_obj )); 
	}
;

recursiveChemInfo:	APDollarSign APOpenParenthesis chemInfo APCloseParenthesis { $$ = new Hold<chem::Root_O>($3->_obj); }
	    ;

intNumber: APNumber { $$ = atoi($1); };




    %%

    //
    //
    //	Epilogue
    //
    //


    typedef	struct {
    char	rchar;  
    char* name ; 
    int	token ; 
    }	SpecialCharacterEntry;

#if 0 
    SpecialCharacterEntry specialCharacters[] = {
    { ';', "APOperatorAndLow", APOperatorAndLow },
    { ',', "APOperatorOr", APOperatorOr },
    { '!', "APOperatorNot", APOperatorNot },
    { '&', "APOPeratorAndHigh", APOperatorAndHigh },
    { '+', "APPositiveCharge", APPositiveCharge },
    //{ '*', "APWildcardOrAtomicMass", APWildcardOrAtomicMass },   
    { '*', "APWildCard", APWildCard},   
    { '=', "APDoubleBond", APDoubleBond },
    { ':', "APAromaticBond", APAromaticBond },
    { '~', "APAnyBond", APAnyBond },
    { '$', "APDollarSign", APDollarSign},
    { '(', "APOpenParenthesis", APOpenParenthesis },
    { ')', "APCloseParenthesis", APCloseParenthesis },
    { '[', "APOpenBracket", APOpenBracket },
{ ']', "APCloseBracket", APCloseBracket },
{ 'D', "APDegree", APDegree },
{ 'H', "APTotalHCount", APTotalHCount },
{ 'h', "APImplicitHCount", APImplicitHCount },
{ 'i', "APPiBondOrbital", APPiBondOrbital },
{ 'Q', "APHeavyAtomTotalBond", APHeavyAtomTotalBond },
{ 'R', "APRingMembershipCount", APRingMembershipCount }, 
{ 'r', "APRingSize", APRingSize },
{ 'v', "APValence", APValence },
{ 'X', "APConnectivity", APConnectivity },
{ '\0', "", 0 }
};
#endif


const char*	elementTable[] = 
{		//First the two character Inorganic elements
"Al",
"Ca",		
"Co",
"Cu",
"Fe",
""		// Last entry
};


const char*	organicElementTable[] = 
{
"Br",       	//Put 2 chars elements on top, or element won't match correctly 
"Cl",
"B",
"b",
"C",
"c",
"F",
"I",
"N",
"n",
"O",
"o",
"P",
"p",
"S",
"s",
""
};

int	msmarts_lex(YYSTYPE* yylval, YYLTYPE* yylloc, msmarts_SParserParam* data)
{
char	ichar;
    do 
    {
	ichar = data->lexStream->get();	 // Read one character
	if (ichar=='\n') 
	return APEndOfLine;
    } while (!data->lexStream->eof() && ichar < ' ' );
    if ( data->lexStream->eof() ) return 0;
#ifdef	DEBUG_ON
    switch (ichar)
    {
        case '&':
	    LEXDPRINT(BF("next char: &amp;"));
	    break;
	case '<':
	    LEXDPRINT(BF("next char: &lt;"));
	    break;
	case '>':
	    LEXDPRINT(BF("next char: &gt;"));
	    break;
	default:
            LEXDPRINT(BF("next char: |%c|") % ichar);
	    break;
    }
#endif


		    //
		    // Handle two and one character elements
		    //
    LEXDPRINT(BF("Checking 1 and 2 chars inorganic elements"));
    for ( int i=0; elementTable[i][0]!='\0'; i++ )
    {
	if ( elementTable[i][0] == ichar ) 
	{
	    if (elementTable[i][1] != '\0') 
	    {
		if ( elementTable[i][1] == data->lexStream->peek() ) 
		    {
		    data->lexStream->get();
		    yylval->eval = elementTable[i];
		    LEXDPRINT(BF("Matched the 2 char inorganic element: |%s|")%elementTable[i]);
		    return APInorganicElement;
		    }
	    } 
	    else 
	    {

		yylval->eval = elementTable[i];
//		    strcpy( yylval->eval , elementTable[i]);
		LEXDPRINT(BF("Matched the 1 char inorganic element: |%s|")%elementTable[i]);	
		return APInorganicElement;
	    }
       }
   }            



// 	Process Both Aliphatic and Aromatic Organic ElementS	
    LEXDPRINT(BF("Checking 1 and 2 chars organic elements"));	
    for (int i=0; organicElementTable[i][0] != '\0' ; i++)
    {	
	if (organicElementTable[i][0] == ichar )   
	{
	    LEXDPRINT(BF("Matched the 1st char in organicElementTable[]"));
	    if (organicElementTable[i][1] != '\0')
	    {
		if ( organicElementTable[i][1] == data->lexStream->peek() )
		{
		    cout<<"Peek char: "<<data->lexStream->peek()<<endl;
		    data->lexStream->get();
//			 strcpy(yylval->eval , organicElementTable[i] );
		    yylval->eval = organicElementTable[i];
		    LEXDPRINT(BF("Matched the 2 char organic element: |%s|")%elementTable[i]);
		    return APOrganicElement;
		}
	    }	
	    else 
	    {
		yylval->eval = organicElementTable[i];
		LEXDPRINT(BF("Matched the 1 char organic element: |%s|\n")
						     % organicElementTable[i]);
						     // If the first character of the element is lowercase then
						     // its an aromatic organic element
		if ( organicElementTable[i][0] >= 'a' && organicElementTable[i][0]<='z') return APAromaticOrganicElement;
		return APOrganicElement;
	    }	
	}
    } 





		    //	
		    // Handle special characters that are one
		    // character long
		    //

    LEXDPRINT(BF("Checking special characters"));	
    switch (ichar) 
    {
    case '*':
	LEXPRINT("APWildCard");
	return APWildCard;
	break;
    case '$':
	LEXPRINT("APDollarSign");
	return APDollarSign;
	break;
    case '(':
	LEXPRINT("APOpenParenthesis");
	return APOpenParenthesis;
	break;
    case ')':
	LEXPRINT("APCloseParenthesis");
	return APCloseParenthesis;
	break;
    case '[':
	LEXPRINT("APOpenBracket");
	return APOpenBracket;
	break;
    case ']':
	LEXPRINT("APCloseBracket");
	return APCloseBracket;
	break;
    case ';':
	LEXPRINT("APOperatorAndLow");
	return APOperatorAndLow;
	break;
    case '?' :
	LEXPRINT("APRingTest");
	return APRingTest;
	break;
    case 'U':	
	// Test if residue of this atom
	// is the same as a previously tagged atom
	LEXPRINT("APResidueTest");
	LEXDPRINT(BF("Identified as APResidueTest"));
	return APResidueTest;
	break;
    case ',':
	LEXPRINT("APOperatorOr");
	return APOperatorOr;
	break;
    case '!':
	LEXPRINT("APOperatorNot");
	return APOperatorNot;
	break;
    case '&':
	LEXPRINT("APOpeartorAndHigh");
	return APOperatorAndHigh;
	break;
    case '~':
	LEXPRINT("APAnyBond");   
	yylval->benum = chem::SABAnyBond;
	return APBond;
	break;	
    case ':' :
	LEXPRINT("APAromaticBond");
	yylval->benum = chem::SABAromaticBond;
	return APBond;
	break;
    case '#':
	if (data->lexStream->peek() >= '0' && data->lexStream->peek() <= '9' ) 
	{
	    LEXPRINT("APAtomicNumber");
	    return APAtomicNumber;
	    break;
	} else if (data->lexStream->peek() == 'G')
	{
	    data->lexStream->get();
	    if (data->lexStream->peek() >= '0' && data->lexStream->peek() <= '9')
	    {
		LEXPRINT("APGroupNumber");
		return APGroupNumber;
		break;
	    }
	    return APError;
	} else if (data->lexStream->peek() == 'T')
	{
	    data->lexStream->get();
	    if (data->lexStream->peek() >= '0' && data->lexStream->peek() <= '9')
	    {
		LEXPRINT("APTransitionMetal");
		return APTransitionMetal;
		break;
	    }
	    return APError;
	} else if (data->lexStream->peek() == 'N')
	{
	    data->lexStream->get();
	    if (data->lexStream->peek() >= '0' && data->lexStream->peek() <= '9')
	    {
		LEXPRINT("APElectronegativeElement");
		return APElectronegativeElement;
		break;
	    }
	    return APError;
	} 
	LEXPRINT("APTripleBond");
	yylval->benum = chem::SABTripleBond;
	return APBond;
	break;
    case '/':
	if (data->lexStream->peek() == '?') 
	{
	    data->lexStream->get();
	    LEXPRINT("APSingleUpOrUnspecified");
	    yylval->benum = chem::SABDirectionalSingleUpOrUnspecified;
	    return  APBond;
	    break;
	}
	LEXPRINT("APDirectionalSingleUp");
	yylval->benum = chem::SABDirectionalSingleUp;
	return  APBond;
	break;
    case '+':
        LEXPRINT("APPositiveCharge");
	return APPositiveCharge;
	break;
    case '-':
        {
	    if (data->lexStream->peek() >= '0' && data->lexStream->peek() <= '9'  ) 
	    {	
		LEXPRINT("APNegativeCharge");
		return APNegativeCharge;
	    }
	    char p = data->lexStream->peek();
	    // if the next char is a ], another - or a logical operation then 
	    // this char represents a negative charge
	    if ( p=='-' || p==']' || p=='&' || p==',' || p==';' )
	    {
		LEXPRINT("APNegativeCharge");
		return APNegativeCharge;
	    }
	    LEXPRINT("APSingleBond");
	    yylval->benum = chem::SABSingleBond;
	    return APBond;
	}
	break;
    case '=':
	LEXPRINT("APDoubleBond");
	yylval->benum = chem::SABDoubleBond;
	return APBond;
	break;	
    case '\\':
	if (data->lexStream->peek() == '?') 
	{
	    data->lexStream->get();
	    LEXPRINT("APSingleDownOrUnspecified");
	    yylval->benum = chem::SABDirectionalSingleDownOrUnspecified;
	    return  APBond;
	    break; 
	}
	LEXPRINT("APDirectionalSingleDown");
	yylval->benum = chem::SABDirectionalSingleDown;
	return  APBond;
	break;
    case '@':
	if ( data->lexStream->peek() == '@' ) {
	    data->lexStream->get(); // pull the second @ out of the stream
	    LEXPRINT("APChiralityClockwise");
	    LEXDPRINT(BF("Identified as APChiralityClockwise"));	
	    return  APChiralityClockwise;
	    break;
	}
	LEXPRINT("APChiralityAntiClockwiseOrAnyRingBond");
	LEXDPRINT(BF("Identified as APChiralityAntiClockwise"));
	return  APChiralityAntiClockwise;
	break;
    case 'D':
	if ( data->lexStream->peek() >= '0' && data->lexStream->peek() <= '9') {
	    LEXPRINT("APDegree");
	    LEXDPRINT(BF("Identified as APDegree"));
	    return APDegree;
	    break;
	}
	return APError;
	break;
    case 'A':
	LEXPRINT("APAliphatic");
	LEXDPRINT(BF("Identified as APAliphatic"));
	return APAliphatic;
	break; 
    case 'a':
	LEXPRINT("APAromatic");
	LEXDPRINT(BF("Identified as APAromatic"));
	return APAromatic;
	break; 
    case 'H':
	LEXPRINT("APTotalHCount");
	LEXDPRINT(BF("Identified as APTotalHCount"));
	return  APTotalHCount;
	break; 
    case 'h':
	LEXPRINT("APImplicitHCount");
	LEXDPRINT(BF("Identifed as APImplicitHCount"));
	return  APImplicitHCount;
	break;
    case 'i':
	if (data->lexStream->peek() >= '0' && data->lexStream->peek() <= '9')
	{
	    LEXPRINT("APAromaticPiElectron");
	    LEXDPRINT(BF("Identified as APAromaticPiElectron"));
	    return APAromaticPiElectron;
	    break;
	}
	LEXPRINT("APPiBondOrbital");
	LEXDPRINT(BF("Identified as APPiBondOrbital"));
	return APPiBondOrbital;
	break;
    case 'Q':
	LEXPRINT("APHeavyAtomTotalBond");
	LEXDPRINT(BF("Identified as APHeavyAtomTotalBond"));
	return APHeavyAtomTotalBond;
	break;  
    case 'R':
	LEXPRINT("APRingMembershipCount");
	LEXDPRINT(BF("Identified as APRingMembership"));
	return APRingMembershipCount;
	break;
    case 'r':
	LEXPRINT("APRingSize");
	LEXDPRINT(BF("Identified as APRingSize"));
	return APRingSize;
	break;
    case '<':
        {
	    LEXPRINT("APLambda");
	    //Create a string for the lambda name
	    int charCount = 0;
	    while ( (data->lexStream->peek()>='0' && data->lexStream->peek()<='9')
		    || (data->lexStream->peek()>='A' && data->lexStream->peek()<='Z')
		    || (data->lexStream->peek()>='a' && data->lexStream->peek()<='z') ) 
	    {
	        ichar = data->lexStream->get();
		if ( charCount >= MAX_CVAL )
		{
	            SIMPLE_ERROR(BF("Exceeded max number of digits in Lambda name"));
		}
		yylval->cval[charCount] = ichar;
		charCount++;
	    }
	    ichar = data->lexStream->get();
	    if ( ichar != '>' )
	    {
	        SIMPLE_ERROR(BF("Illegal lambda name must match <[A-Za-z0-9]+> - instead saw[%c]") % ichar );
	    }	
	    yylval->cval[charCount] = '\0';
	    lisp_LOG(BF("APLambda<%s>") % (yylval->cval) );
	    return  APLambda;
	}
	break;
    case 'v':
	LEXPRINT("APValence");
	LEXDPRINT(BF("Identified as APValence"));	
	return APValence;
	break;
    case 'X':
	LEXPRINT("APConnectivity");  
	LEXDPRINT(BF("Identified as APConnectivity"));	
	return APConnectivity;
	break;
    case 'x': // AM1-BCC "x" type
	      // matches to atoms that can be in aromatic rings
	      // C(x3),N(x2),P(x2),N+(x3),P+(x3),O+(x2),S+(x2)
	      // See Jakalian, Jack, Bayly J Comput Chem 23: 1623–1641, 2002
	LEXPRINT("AM1-BCC-x");
	LEXDPRINT(BF("Identified as AM1-BCC-x"));
	return AM1_BCC_x;
	break;
    case 'y': // AM1-BCC "y" type
	      // matches to atoms that can be in 5-member aromatic rings
	      // C-(x2),N-(x2),O(x2),S(x2),N(x3),P(x3)
	      // See Jakalian, Jack, Bayly J Comput Chem 23: 1623–1641, 2002
	LEXPRINT("AM1-BCC-y");
	LEXDPRINT(BF("Identified as AM1-BCC-y"));
	return AM1_BCC_y;
	break;
    }
		    //	
		    // Handle special characters that are two 
		    // characters long
		    //
    if (ichar == 'L' && data->lexStream->peek() == 'p')
    {
        data->lexStream->get();	
        LEXPRINT("APLonePair");
        return APLonePair;
    }

			    // Parse a number
    if ( ichar>='0' && ichar<='9' ) 
    {
	    //Create a string for number	
	    int digitCount = 0;
	    yylval->cval[digitCount] = ichar;
	    digitCount++;
	    while ( data->lexStream->peek()>='0' && data->lexStream->peek()<='9' ) 
	    {
		    ichar = data->lexStream->get();
		    if ( digitCount >= MAX_CVAL )
		    {
		        SIMPLE_ERROR(BF("Exceeded max number of digits in Number"));
		    }
		    yylval->cval[digitCount] = ichar;
		    digitCount++;
	    }
	    yylval->cval[digitCount] = '\0';
	    lisp_LOG(BF("APNumber = (%s)") % (yylval->cval) );
	    return  APNumber;
    }
    

    LEXPRINT("APError");
    return APError;


}



#if 0
int main()
{
    cout<<"Input CHEM_INFO : ";
    data->lexStream = &cin;

    yyparse();


return 0;
}
#endif


chem::SmartsRoot_sp smarts_compile(const string& input, core::Lisp_sp lisp, stringstream& errorStream)
{_G();
    msmarts_SParserParam p;
    stringstream sin(input);
    p.lisp = lisp;
    p.expression = _Nil<SmartsRoot_O>();
    p.lexStream = &sin;
    p.msmartsErrorStream = &errorStream;
    lisp_LOG(BF("Entering msmarts_parse"));
    if ( msmarts_parse(&p) )
    {
	SIMPLE_ERROR(BF("%s") % (p.msmartsErrorStream->str()));
    }
    return p.expression;
}


void msmarts_error(YYLTYPE* yyloc, msmarts_SParserParam* data, const char* message )
{
    *(data->msmartsErrorStream) << "Error: " << message << endl;
}
    


