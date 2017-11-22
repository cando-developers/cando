       
       
       
//
//
//   	antechamberParse.yy  
//
//
%{
#define DEBUG_LEVEL_FULL
#include <string.h>
#include <iostream>
#include <string>
#include <vector>
#include <istream>

#include <clasp/core/foundation.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/chem/hold.h>

#include <cando/chem/chemInfo.h>
#include <cando/chem/gaff_TypeParser.h>
#include <cando/chem/gaff_ParserParam.h>
#include <cando/chem/gaff_LocationType.h>


using namespace std;
using namespace core;



#define	LEXDEBUG	1

#ifdef	LEXDEBUG
#define	LEXPRINT(x) {LOG(BF("Token: %s") % x );}
#define	LEXDPRINT(x) { LOG(BF("%s") % (x));}
#else
#define	LEXPRINT(x)
#define	LEXDPRINT(x)
#endif

#define	_BT()	_BLOCK_TRACE("block trace")

#define	GP() (data)

void gaff_error(YYLTYPE* yyloc, gaff_SParserParam* data, const char* message );

int	gaff_lex(YYSTYPE* yylval, YYLTYPE* yylloc, gaff_SParserParam* data);



%}

// %define api.prefix {gaff_}
// %n a m e - p r e f i x ="gaff_"
%define api.pure full
%locations
%lex-param {gaff_SParserParam* data}
%parse-param {gaff_SParserParam* data}
%define parse.trace





%token APEndOfLine
%token <carr> APElement
%token <carr> APElementWildCard
%token APWildCard;
%token APOpenParenthesis
%token APCloseParenthesis
%token APOpenBracket
%token APOpenAngleBracket
%token APCloseBracket
%token APCloseAngleBracket
%token APOperatorAnd
%token APOperatorOr
%token APColon
%token <ival> APAntechamberARLevel
%token APAntechamberRingMembership
%token APAntechamberNoRingMembership
%token <benum> APBond
%token <ival> APNumber
%token <carr> APTag
%token APError
%token APNotBondedToPrevious
%token APBondedToPrevious
%token  APAtomTypeDefinition

%left	APOperatorAnd
%left	APOperatorOr
%left	APOperatorNot



%type <carr>    tagNameOrNull
%type <carr>    element;
%type <carr>    tagOrElement;
%type <carr>    typeName;
%type <ival>    focusAtomicNumber;
%type <ival>    focusNumberOfAttachedAtoms;
%type <ival>    focusNumberOfAttachedHydrogenAtoms;
%type <ival>    focusNumberOfElectronWithdrawingAtoms;
//%type <antechamberFocusAtomMatch>    focusTest;
//%type <antechamberRoot>	antechamberMatch;
%type <antechamberRoot>	antechamberLine;
%type <bondListMatchNode>    chemicalEnvironment;
%type <afterMatchBondTest>    bondDefinitions ;
%type <residueList>    residueList;
%type <logical>    bracketedAtomicProperty
%type <logical>    logOp 
%type <atomTest>    atomicTest
%type <bondListMatchNode>    branches
%type <bondListMatchNode>    branchesList
%type <bondListMatchNode>    chain 
%type <antechamberBondTest>    chemicalEnvironmentAtom
%type <logical>    bracketedAtomPropertyOrNull

//%destructor { if ($$) chem::deleteHolder($$); } 	antechamberMatch;
%destructor { if ($$) chem::deleteHolder($$); } 	antechamberLine;
%destructor { if ($$) chem::deleteHolder($$); }     chemicalEnvironment;
%destructor { if ($$) chem::deleteHolder($$); }     bondDefinitions ;
%destructor { if ($$) chem::deleteHolder($$); }     residueList;
%destructor { if ($$) chem::deleteHolder($$); }     bracketedAtomicProperty
%destructor { if ($$) chem::deleteHolder($$); }     logOp 
%destructor { if ($$) chem::deleteHolder($$); }     atomicTest
%destructor { if ($$) chem::deleteHolder($$); }     branches
%destructor { if ($$) chem::deleteHolder($$); }     branchesList
%destructor { if ($$) chem::deleteHolder($$); }     chain 
%destructor { if ($$) chem::deleteHolder($$); }     chemicalEnvironmentAtom
%destructor { if ($$) chem::deleteHolder($$); }     bracketedAtomPropertyOrNull

%%




antechamberMatch: 
        APAtomTypeDefinition 
                antechamberLine { data->expression = $2->value(); LOG(BF("Got AntechamberMatch tree=%0X") % $2 ); }
        ;

antechamberLine: 
          typeName
                residueList 
                focusAtomicNumber 
                focusNumberOfAttachedAtoms 
                focusNumberOfAttachedHydrogenAtoms 
                focusNumberOfElectronWithdrawingAtoms 
                bracketedAtomicProperty 
                chemicalEnvironment 
                bondDefinitions 
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( $2->value(), $3, $4, $5, $6, $7->value() );
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus,$8->value(),$9->value())); 
                }
        | typeName
                residueList 
                focusAtomicNumber 
                focusNumberOfAttachedAtoms 
                focusNumberOfAttachedHydrogenAtoms 
                focusNumberOfElectronWithdrawingAtoms 
                bracketedAtomicProperty 
                chemicalEnvironment 
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( $2->value(), $3, $4, $5, $6, $7->value() );
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus,$8->value())); 
                }
        | typeName
                residueList 
                focusAtomicNumber 
                focusNumberOfAttachedAtoms 
                focusNumberOfAttachedHydrogenAtoms 
                focusNumberOfElectronWithdrawingAtoms 
                bracketedAtomicProperty 
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( $2->value(), $3, $4, $5, $6, $7->value() );
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus)); 
                }
        | typeName
                residueList 
                focusAtomicNumber 
                focusNumberOfAttachedAtoms 
                focusNumberOfAttachedHydrogenAtoms 
                focusNumberOfElectronWithdrawingAtoms 
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    core::T_sp temp = _Nil<core::T_O>();
                    focus = chem::AntechamberFocusAtomMatch_O::create( $2->value(), $3, $4, $5, $6, temp );
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus)); 
                }
        | typeName
                residueList 
                focusAtomicNumber 
                focusNumberOfAttachedAtoms 
                focusNumberOfAttachedHydrogenAtoms 
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( $2->value(), $3, $4, $5, -1);
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus)); 
                }
        | typeName
                residueList 
                focusAtomicNumber 
                focusNumberOfAttachedAtoms 
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
		    focus = chem::AntechamberFocusAtomMatch_O::create( $2->value(), $3, $4, -1, -1);
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus));
                }
        | typeName
                residueList 
                focusAtomicNumber 
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( $2->value(), $3, -1, -1, -1 );
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus)); 
                }
        | typeName
                residueList 
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( $2->value(), -1, -1, -1, -1);
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus)); 
                }
        | typeName
                APEndOfLine
                {
		    _BLOCK_TRACE("Got antechamber line");
                    gc::Nilable<chem::AtomOrBondMatchNode_sp> focus = _Nil<core::T_O>();
                    $$ = chem::createHolder<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create(chem::chemkw_intern(std::string($1)),focus));
                }
        ;

typeName:
         tagOrElement { strcpy($$,$1); }
        | tagOrElement APOperatorOr APNumber { stringstream ss; ss << $1 << "." << $3; strcpy($$,ss.str().c_str()); }
        | tagOrElement APOperatorOr APNumber tagOrElement { stringstream ss; ss << $1 << "." << $3 << $4; strcpy($$,ss.str().c_str()); }
        | tagOrElement APOperatorOr tagOrElement { stringstream ss; ss << $1 << "." << $3; strcpy($$,ss.str().c_str()); }
        ;

tagOrElement:
          APTag   { strcpy($$,$1); }
        | APElement { strcpy($$,$1); }
        ;





residueList:    APWildCard      { $$ = chem::createHolder<chem::ResidueList_O>(chem::ResidueList_O::create()); }
        ;
       // Later add residue list parsing

focusAtomicNumber: 
          APWildCard   { $$ = -1; }
        | APNumber { $$ = $1; }
        ;

focusNumberOfAttachedAtoms : 
          APWildCard { $$ = -1; }
        | APNumber { $$ = $1; }
        ;

focusNumberOfAttachedHydrogenAtoms:
          APWildCard { $$ = -1; }
        | APNumber { $$ = $1; }
        ;

focusNumberOfElectronWithdrawingAtoms: 
          APWildCard { $$ = -1; }
        | APNumber { $$ = $1; }
        ;

bracketedAtomicProperty: 
          APWildCard     { _BT(); $$ = chem::createHolder<chem::Logical_O>(); }
        | APOpenBracket logOp APCloseBracket { _BT(); 
		LOG(BF("Got bracketedAtomicProperty") ); // vp0(( "Got bracketedAtomicProperty" ));
		$$ = chem::createHolder<chem::Logical_O>($2->value()); }
        ;


chemicalEnvironment: 
          APWildCard      { _BT(); $$ = chem::createHolder<chem::BondListMatchNode_O>(); }
        | branches { _BT(); $$ = chem::createHolder<chem::BondListMatchNode_O>($1->value()); }
        ;



bondDefinitions:        
          APTag APColon APTag APColon APBond {_BT(); $$ = chem::createHolder<chem::AfterMatchBondTest_O>(chem::AfterMatchBondTest_O::create( chem::chemkw_intern(std::string($1)), chem::chemkw_intern(std::string($3)), $5 )); }
        ;


logOp:  atomicTest {_BT(); $$ = chem::createHolder<chem::Logical_O>(chem::Logical_O::create( chem::logIdentity, $1->value() )); }
        | logOp APOperatorAnd logOp { _BT();$$ = chem::createHolder<chem::Logical_O>(chem::Logical_O::create(chem::logHighPrecedenceAnd, $1->value(), $3->value() )); }
        | logOp APOperatorOr logOp {_BT(); $$ = chem::createHolder<chem::Logical_O>(chem::Logical_O::create(chem::logOr, $1->value(), $3->value() )); }
        ;

atomicTest:     APBond  { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPInBond, $1, 1)); }
        | APNumber APBond  { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPInBond, $2, $1)); }
        | APBond APNotBondedToPrevious { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPNotBondedToPrevious,$1)); }
        | APBond APBondedToPrevious { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPBondedToPrevious,$1)); }
        | APAntechamberARLevel { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPArLevel, $1 )); }
        | APAntechamberRingMembership { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPRingSize, 1 )); }
        | APAntechamberRingMembership APNumber { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPRingSize, $2 )); }
        | APNumber APAntechamberRingMembership APNumber { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPRingSize, $3, $1 )); }
        | APAntechamberNoRingMembership { _BT(); $$ = chem::createHolder<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPNoRing)); }
        ;


branches: APOpenParenthesis branchesList APCloseParenthesis { _BT(); $$=$2; }
        ;

branchesList:    chain { _BT(); $$ = chem::createHolder<chem::BondListMatchNode_O>($1->value()); }
        | branchesList APOperatorAnd chain { _BT(); $$ = chem::createHolder<chem::BondListMatchNode_O>(chem::Branch_O::create($1->value(),$3->value())); }
        ;


chain: chemicalEnvironmentAtom { _BT();$$ = chem::createHolder<chem::BondListMatchNode_O>(chem::Chain_O::create($1->value()));}
        | chemicalEnvironmentAtom branches { _BT(); $$ = chem::createHolder<chem::BondListMatchNode_O>(chem::Chain_O::create($1->value(),$2->value())); }
        ;


chemicalEnvironmentAtom: 
        element 
                APNumber 
                bracketedAtomPropertyOrNull 
                tagNameOrNull { _BT(); 
                    $$ = chem::createHolder<chem::AntechamberBondTest_O>(chem::AntechamberBondTest_O::create(chem::chemkw_intern(std::string($1)),$2,$3->value(),chem::chemkw_intern(std::string($4)))); 
                }
        | element 
                bracketedAtomPropertyOrNull 
                tagNameOrNull { _BT(); 
                    $$ = chem::createHolder<chem::AntechamberBondTest_O>(chem::AntechamberBondTest_O::create(chem::chemkw_intern(std::string($1)), -1, $2->value(), chem::chemkw_intern(std::string($3))));
                }
        ;

element:
          APElement { _BT(); strcpy($$,$1); }
        | APElementWildCard { _BT(); strcpy($$,$1); }
        ;

bracketedAtomPropertyOrNull: /* empty */        { _BT(); $$ = chem::createHolder<chem::Logical_O>(); }
        | APOpenBracket logOp APCloseBracket { _BT(); $$ = chem::createHolder<chem::Logical_O>($2->value()); }
        ;


tagNameOrNull :       /* empty */       { _BT(); strcpy($$,""); }
        | APOpenAngleBracket APTag APCloseAngleBracket { _BT(); strcpy($$,$2); }
        ;





%%

//
//
//	Epilogue
//
//



const char*	acElementTable[] =  {
"He",
"Li",
"Be",
"Ne",
"Na",
"Mg",
"Al",
"Si",
"Cl",
"Ar",
"Ca",
"Sc",
"Ti",
"Cr",
"Mn",
"Fe",
"Co",
"Ni",
"Cu",
"Zn",
"Ga",
"Ge",
"As",
"Se",
"Br",
"Kr",
"Rb",
"Sr",
"Zr",
"Nb",
"Mo",
"Tc",
"Ru",
"Rh",
"Pd",
"Ag",
"Cd",
"In",
"Sn",
"Sb",
"Te",
"Xe",
"Cs",
"Ba",
"La",
"Ce",
"Pr",
"Nd",
"Pm",
"Sm",
"Eu",
"Gd",
"Tb",
"Dy",
"Ho",
"Er",
"Tm",
"Yb",
"Lu",
"Hf",
"Ta",
"Re",
"Os",
"Ir",
"Pt",
"Au",
"Hg",
"Tl",
"Pb",
"Bi",
"Po",
"At",
"Rn",
"Fr",
"Ra",
"Ac",
"Th",
"Pa",
"Np",
"Pu",
"Am",
"Cm",
"Bk",
"Cf",
"Es",
"Fm",
"Md",
"No",
"Lr",
"H",
"W",
"U",
"B",
"C",
"N",
"O",
"F",
"P",
"S",
"K",
"V",
"I",
"Y",
"" };



void gaff_error(YYLTYPE* yyloc, gaff_SParserParam* data, const char* message )
{
    *(data->gaffErrorStream) << "Error: " << message << endl;
}


int	gaff_lex(YYSTYPE* gaff_lval, YYLTYPE* gaff_lloc, gaff_SParserParam* data )
{ 
char	ichar;
    do {
	ichar = data->acLexStream->get();	 // Read one character
	if (ichar=='\n') 
	return APEndOfLine;
    } while (!data->acLexStream->eof() && ichar <= ' ' );
    if ( data->acLexStream->eof() ) return 0;

                //
                // Check wildcards
                //
    if ( data->wildElementDictionary.notnilp() ) {
        if ( data->wildElementDictionary->lexWildNameMatches2Char(ichar,
                                        data->acLexStream->peek()) )
	{
            int ichar2 = data->acLexStream->get();
            gaff_lval->carr[0] = ichar;
            gaff_lval->carr[1] = ichar2;
            gaff_lval->carr[2] = '\0';
	    LEXDPRINT(("Matched the 2 char element wildcard: |%s|",
                                gaff_lval->carr));
            return APElementWildCard;
        }
        if ( data->wildElementDictionary->lexWildNameMatches1Char(ichar) ) {
            gaff_lval->carr[0] = ichar;
            gaff_lval->carr[1] = '\0';
	    LEXDPRINT(("Matched the 1 char element wildcard: |%s|",
                                gaff_lval->carr));
            return APElementWildCard;
        }
    }
        

		    //	
		    // Handle special characters that are one
		    // character long
		    //

//    LEXDPRINT(("Checking special characters for ichar(%c)", ichar));	
    switch (ichar) {
	    case '*':
		    LEXPRINT("APWildCard");
		    return APWildCard;
		    break;
	    case '&':
		    LEXPRINT("APEndOfLine");
		    return APEndOfLine;
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
	    case '<':
		    LEXPRINT("APOpenAngleBracket");
		    return APOpenAngleBracket;
		    break;
	    case '>':
		    LEXPRINT("APCloseAngleBracket");
		    return APCloseAngleBracket;
		    break;
	    case ',':
		    LEXPRINT("APOperatorAnd");
		    return APOperatorAnd;
		    break;
	    case '.':
		    LEXPRINT("APOperatorOr");
		    return APOperatorOr;
		    break;
            case ':':
		    LEXPRINT("APColon");
		    return APColon;
		    break;
    }

		    //	
		    // Handle special characters that are two or three 
		    // characters long
		    //
    LOG(BF("peek(%c)") % (acLexStream->peek() ) ); // vp0(( "peek(%c)", acLexStream->peek() ));
    LOG(BF("Checking ATD") ); // vp0(("Checking ATD"));
    if ( ichar == 'A' && data->acLexStream->peek() == 'T' ) {
        data->acLexStream->get();
        LOG(BF("peek2(%c)") % (acLexStream->peek() ) ); // vp0(( "peek2(%c)", acLexStream->peek() ));
        if ( data->acLexStream->peek() != 'D' ) {
            return APError;
        }
        data->acLexStream->get();
        LEXPRINT("APAtomTypeDefinition");
        return APAtomTypeDefinition;
    }
        
    LOG(BF("Checking AR") ); // vp0(("Checking AR"));
    if (ichar == 'A' && data->acLexStream->peek() == 'R') {
        data->acLexStream->get();
        int aridx = data->acLexStream->get()-'0';
        if ( aridx >=1 && aridx <=5 ) {
            gaff_lval->ival = aridx;
	    LOG(BF("AR level=%d") % (aridx ) ); // vp0(( "AR level=%d", aridx ));
            return APAntechamberARLevel;
        }
        LEXPRINT("APAntechamberARLevel");
        return APError;
    }

    LOG(BF("Checking RG") ); // vp0(("Checking RG"));
    if (ichar == 'R' && data->acLexStream->peek() == 'G') {
        data->acLexStream->get();	
        LEXPRINT("APAntechamberRingMembership");
        return APAntechamberRingMembership;
    }

    LOG(BF("Checking NG") ); // vp0(("Checking NG"));
    if (ichar == 'N' && data->acLexStream->peek() == 'G') {
        data->acLexStream->get();	
        LEXPRINT("APAntechamberNoRingMembership");
        return APAntechamberNoRingMembership;
    }

    LOG(BF("Checking SB") ); // vp0(("Checking SB"));
    if (ichar == 'S' && data->acLexStream->peek() == 'B') {
        data->acLexStream->get();	
        LEXPRINT("APSingleBond");
        gaff_lval->benum = chem::SABSingleBond;
        return APBond;
    }
    if (ichar == 'D' && data->acLexStream->peek() == 'B') {
        data->acLexStream->get();	
        LEXPRINT("APDoubleBond");
        gaff_lval->benum = chem::SABDoubleBond;
        return APBond;
    }
    if (ichar == 'D' && data->acLexStream->peek() == 'L') {
        data->acLexStream->get();	
        LEXPRINT("APDelocalizedBond");
        gaff_lval->benum = chem::SABDelocalizedBond;
        return APBond;
    }
    if (ichar == 'T' && data->acLexStream->peek() == 'B') {
        data->acLexStream->get();	
        LEXPRINT("APTripleBond");
        gaff_lval->benum = chem::SABTripleBond;
        return APBond;
    }
    if (ichar == 's' && data->acLexStream->peek() == 'b') {
        data->acLexStream->get();	
        LEXPRINT("APSingleOrAromaticBond");
        gaff_lval->benum = chem::SABSingleOrAromaticBond;
        return APBond;
    }
    if (ichar == 'd' && data->acLexStream->peek() == 'b') {
        data->acLexStream->get();	
        LEXPRINT("APDoubleOrAromaticBond");
        gaff_lval->benum = chem::SABDoubleOrAromaticBond;
        return APBond;
    }
    if (ichar == 't' && data->acLexStream->peek() == 'b') {
        data->acLexStream->get();	
        LEXPRINT("APTripleBond");
        gaff_lval->benum = chem::SABTripleBond;
        return APBond;
    }
    if (ichar == 'a' && data->acLexStream->peek() == 'n') {
        data->acLexStream->get(); // eat 'n'
        char fc = data->acLexStream->get(); // eat 'y'
        if ( fc != 'y' ) {
            return APError;
        }
        LEXPRINT("APAnyBond");
        gaff_lval->benum = chem::SABAnyBond;
        return APBond;
    }

    if ( ichar == '\'' ) {
        if ( data->acLexStream->peek() == '\'' ) {
            data->acLexStream->get();
            return APNotBondedToPrevious;
        }
        return APBondedToPrevious;
    }
		    //
		    // Handle two and one character elements
		    //
    LEXDPRINT(("Checking 2 and 1 chars elements"));	
    for ( int i=0; acElementTable[i][0]!='\0'; i++ ) {
	if ( acElementTable[i][0] == ichar ) {
	    if (acElementTable[i][1] != '\0') {
		if ( acElementTable[i][1] == data->acLexStream->peek() ) {
		    data->acLexStream->get();
		    LEXDPRINT(("Matched the 2 char element: |%s|",
                                                acElementTable[i]));	
		    strcpy( gaff_lval->carr ,  acElementTable[i] );
//		    gaff_lval->eval = acElementTable[i];
		    return APElement;
	        }
	    } else {
		LEXDPRINT(("Matched the 1 char element:|%s|",
                                                acElementTable[i]));
                strcpy( gaff_lval->carr, acElementTable[i] );
//		gaff_lval->eval = acElementTable[i];
		return APElement;
	    }
       }
   }

    if ((ichar>='a'&&ichar<='z')||(ichar>='A'&&ichar<='Z')) {
        string tagStr;
        tagStr = ichar;
        while ( (data->acLexStream->peek()>='a'&&data->acLexStream->peek()<='z') ||
                (data->acLexStream->peek()>='A'&&data->acLexStream->peek()<='Z') ||
                (data->acLexStream->peek()>='0'&&data->acLexStream->peek()<='9') ) {
            ichar = data->acLexStream->get();
            tagStr += ichar;
        }
        strcpy(gaff_lval->carr,tagStr.c_str());
        LOG(BF("Got tag(%s)") % (gaff_lval->carr) ); // vp0(( "Got tag(%s)",gaff_lval->carr));
        LEXPRINT("APTag");
        return APTag;
    }



			    // Parse a number
    if ( ichar>='0' && ichar<='9' ) 
    {
	    //Create a string for number	
	    string numStr;
	    numStr = ichar;
	    while ( data->acLexStream->peek()>='0' && data->acLexStream->peek()<='9' ) 
	    {
		    ichar = data->acLexStream->get();
		    numStr += ichar;
	    }
	    
	    gaff_lval->ival = atoi(numStr.c_str());
	    LEXDPRINT(("APNumber(%d)", gaff_lval->ival));
	    return  APNumber;
    }
    

    LEXPRINT("APError");
    return APError;


}





chem::AntechamberRoot_sp gaff_compile(const string& input, 
				chem::WildElementDict_sp dict,
				stringstream& errorStream)
{_G();
    gaff_SParserParam p;
    stringstream sin(input);
    p.expression = chem::AntechamberRoot_O::create();
    p.wildElementDictionary = dict;
    p.acLexStream = &sin;
    p.gaffErrorStream = &errorStream;
    LOG(BF("Entering gaffparse"));
    if ( gaff_parse(&p) )
    {
	SIMPLE_ERROR(BF("%s") % (p.gaffErrorStream->str()));
    }
    return p.expression;
}



    

