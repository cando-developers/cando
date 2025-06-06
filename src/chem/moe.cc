/*
    File: moe.cc
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
       
#define DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//





/*
 *	moe.cc
 *
 *	Functions that read and write MOE files.
 *
 */
#include <string.h>
#include <clasp/core/common.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <queue>
#include <stdio.h>
#include <map>
#include <cando/chem/matter.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/chem/loop.h>
#include <cando/chem/moe.h>
#include <cando/chem/restraint.h>


namespace chem {
    typedef	std::queue<string>	QueueString;
};







namespace chem {


class	ParaPrimValue {
public:
	virtual ~ParaPrimValue() {};
};

template <class TYPE>
class	ParaValue: public ParaPrimValue {
public:
	TYPE	val;

		ParaValue(TYPE a) {this->val=a;};
	void	setValue(TYPE v) {this->val = v;};
	TYPE&	getValue()	{return this->val;};
};


class	ParaVector: public ParaPrimValue {
public:
	vector<ParaPrimValue*>	vals;

	virtual ~ParaVector()
	{
//	vector<ParaPrimValue*>::iterator	i;
	    this->erase();
//	    for (i=this->vals.begin();i!=this->vals.end();i++) {
//		destroy (*i);
//	    }
	}
	void	push_back(ParaPrimValue* p) {this->vals.push_back(p);};
	ParaPrimValue*	operator[](int i)	{return this->vals[i];};
	int	size()	{return this->vals.size();};
	void	erase()	{this->vals.erase(this->vals.begin(),this->vals.end());};
};
template <class TYPE>
class	ParaArray: public ParaPrimValue {
public:
	vector<TYPE>	vals;
	void	push_back(TYPE p) {this->vals.push_back(p);};
	TYPE	operator[](int i)	{return this->vals[i];};
	int	size()	{return this->vals.size();};
	void	erase()	{this->vals.erase(this->vals.begin(),
						this->vals.end());};
};




typedef	ParaValue<int>		ParaInteger;
typedef	ParaValue<int>		ParaHexInteger;
typedef	ParaValue<string>	ParaString;
typedef	ParaValue<string>	ParaToken;
typedef	ParaValue<double>	ParaReal;
typedef	ParaValue<string>	ParaChar;
typedef	ParaArray<int>		ParaArrayInt;
typedef	ParaArray<double>	ParaArrayReal;



class	ParaTable {
public:
	map<string,int>		fieldIndexes;
	vector<ParaVector>	vectors;

	void	addFieldVector(string f,ParaVector& v) { this->fieldIndexes[f]=this->vectors.size(); this->vectors.push_back(v); };

	ParaVector&	getVector(int f) {return this->vectors[f];};

	ParaVector&	getFieldVector(string f) {return this->vectors[this->fieldIndexes[f]];};
	int	getFieldIndex(string f) {return this->fieldIndexes[f];};




	void	erase()	{ this->vectors.erase(this->vectors.begin(),
					this->vectors.end());
			  this->fieldIndexes.erase(this->fieldIndexes.begin(),
						this->fieldIndexes.end() );};
};





#define	MAXWORD	1024

class	MoeReadFile {
private:
	string		fileName;
	bool		mustDestroy_fMoe;
    std::istream		*fMoe;
	bool		done;
	char		nextLine[MAXWORD];

		// Paragraph stuff
	string		moeFileVersion;
	string		paraName;
	int		paraRows;
	vector<string>	paraFieldNames;
	vector<string>	paraFieldTypes;
    std::queue<string>	body;

	void		parseHeader();
	void		splitLine( QueueString& q, char* line );
	void		readNextLine();

public:
	~MoeReadFile();
	void		openFileName( const char* sFileName );
    void		openStream(std::istream* fIn);
	bool		nextParaTableExists() {return !this->done;};
	string		nextParaTableName() {return this->paraName;};
	int		nextParaTableRows() {return this->paraRows;};
	int		nextParaTableFieldCount() {return this->paraFieldNames.size();};
	string		nextParaTableFieldName(int i) {return this->paraFieldNames[i];};
	string		nextParaTableFieldType(int i) {return this->paraFieldTypes[i];};
	void		readParaTable( ParaTable* p );
	void		skipParaTable();

};




///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////



//
//	_moeWriteRestraint(fout,restraint,agg);
//
    void	_moeWriteRestraint( std::ostream& fout, Restraint_sp rest,
					Aggregate_sp agg )
{
	// At this point all restraints must be absolute
    if ( rest.isA<RestraintDistance_O>() )
    {
      SIMPLE_ERROR("Add support for distance restraints");
#if 0
      RestraintDistance_sp distRest = (rest).as<RestraintDistance_O>();
	fout << "distance 2 "
		<< distRest->getAtomA()->getTempInt() << " "
		<< distRest->getAtomB()->getTempInt() << " "
		<< " 2 "
		<< distRest->getMin() << " "
		<< distRest->getMax() << " "
             << distRest->getWeight() << std::endl;
	return;
#endif
    }
    if ( rest.isA<RestraintAngle_O>() )
    {
      RestraintAngle_sp distAngle = (rest).as<RestraintAngle_O>();
	fout << "angle 3 "
		<< distAngle->getAtomA()->getTempInt() << " "
		<< distAngle->getAtomB()->getTempInt() << " "
		<< distAngle->getAtomC()->getTempInt() << " "
		<< " 2 "
		<< distAngle->getMin() << " "
		<< distAngle->getMax() << " "
             << distAngle->getWeight() << std::endl;
	return;
    }
    RestraintDihedral_sp distDih = (rest).as<RestraintDihedral_O>();
    fout << "dihedral 4 "
	    << distDih->getAtomA()->getTempInt() << " "
	    << distDih->getAtomB()->getTempInt() << " "
	    << distDih->getAtomC()->getTempInt() << " "
	    << distDih->getAtomD()->getTempInt() << " "
	    << " 2 "
	    << distDih->getDegrees() << " "
         << distDih->getWeight() << std::endl;
    return;
}





//
//	isMoePrefixChar
//
//	Return true if the character is a Moe prefix character
//
bool	isMoePrefixChar( char c )
{
    if ( c == '$' ) return true;
    if ( c == '%' ) return true;
    if ( c == '&' ) return true;
    if ( c == '~' ) return true;
    return false;
}



//
//	isMoeEncodeChar
//
//	Return true if the character is a valid MOE file character
static bool	isMoeEncodeChar( char c )
{
    if ( c >= 'A' && c <= 'Z' ) return true;
    if ( c >= 'a' && c <= 'z' ) return true;
    if ( c >= '0' && c <= '9' ) return true;
    if ( c == '+' ) return true;
    if ( c == '/' ) return true;
    return false;
}

//
//	isMoeRawChar
//
//	Return true if the character is a valid MOE file character
static bool	isMoeRawChar( char c )
{
    if ( isMoeEncodeChar(c) ) return true;
    if ( c == '(' ) return true;
    if ( c == ')' ) return true;
    if ( c == '-' ) return true;
    if ( c == '*' ) return true;
    if ( c == '.' ) return true;
    if ( c == ',' ) return true;
    if ( c == ':' ) return true;
    if ( c == ';' ) return true;
    if ( c == '=' ) return true;
    if ( c == '?' ) return true;
    if ( c == '[' ) return true;
    if ( c == '_' ) return true;
    if ( c == ']' ) return true;
    if ( c == '{' ) return true;
    if ( c == '|' ) return true;
    if ( c == '}' ) return true;
    if ( c == '^' ) return true;
    if ( c == ' ' ) return true;
    return false;
}


//
//	isMoeValidChar
//
//	Return true if the character is a valid MOE file character
static bool	isMoeValidChar( char c )
{
    if ( isMoeRawChar(c) ) return true;
    if ( isMoePrefixChar(c) ) return true;
    if ( c == '!' ) return true;
    return false;
}





//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//
//	MoeReadFile
//
//


MoeReadFile::~MoeReadFile()
{
    if ( this->mustDestroy_fMoe) {
	delete this->fMoe;
    }
}

//
//	splitLine
//
//	Split a line of whitespace separated words into separate
//	strings and put them in a queue.
void	MoeReadFile::splitLine( QueueString& q, char* line)
{
    char*	start;
    char	caWord[MAXWORD];
    char*	copy;
    char	prefix, encode;
    string	retVal;

    start = line;
    while (*start) {
	LOG("Parsing rest of line: |{}|" , start  );
        copy = caWord;
	// skip whitespace
	while ( (*start) && (isspace(*start)) ) start++;
	// copy the characters into caWord
	while ( (*start) && (!isspace(*start)) ) {
	    // Skip invalid characters
	    if ( !isMoeValidChar(*start) ) {
		start++;
		continue;
	    }
	    // Handle encoded characters
	    if ( isMoePrefixChar(*start) ) {
		prefix = *start;
		start++;
		encode = *start;
		if ( isMoeEncodeChar(encode) ) {
		    *copy = ((prefix - 1)<<6) | ((encode-1)&0x3F);
		    copy++;
		} else {
		    *copy = prefix;
		    copy++;
		}
		continue;
	    }
	    // Convert '!' to spaces
	    if ( *start=='!' ) {
		*copy = ' ';
		copy++;
		start++;
		continue;
	    } else {
		*copy = *start;
		copy++;
		start++;
		continue;
	    }
	}
	*copy = '\0';

	// Now we have an uncoded word in caWord
	// convert it into a string
	if ( strcmp(caWord,"$") == 0 ) {
	    retVal = "";
	} else {
	    retVal = caWord;
	}
	LOG("Word: |{}|" , retVal.c_str()  );
	q.push(retVal);
    }
}


//
//	parseHeader
//
//	Read the header and parse the strings into the header queue.
//
void	MoeReadFile::parseHeader()
{
char		line[MAXWORD];
char		buf[MAXWORD];
std::queue<string>	words;
string		word;

    strcpy( line, this->nextLine );
    if ( line[0] != '#' ) {
	strcpy( buf, "Non moe header: " );
	strcat( buf, line );
	SIMPLE_ERROR("Invalid moe header");
    }
    this->splitLine( words, line );
    this->paraName = words.front();
    this->paraRows = 0;
    words.pop();
    if ( this->paraName == "moe" ) {
	this->moeFileVersion = words.front();
	words.pop();
	goto DONE;
    }
		// Get the number of rows
    if ( !words.empty() ) {
	word = words.front();
	words.pop();
	strncpy( buf, &(word[0]), word.length() );
	buf[word.length()] = '\0';
	this->paraRows = atoi(buf);
    }
		// Parse out the field names and types
    this->paraFieldNames.erase(this->paraFieldNames.begin(),
				this->paraFieldNames.end() );
    this->paraFieldTypes.erase(this->paraFieldTypes.begin(),
				this->paraFieldTypes.end() );
    while ( ! words.empty() ) {
	this->paraFieldNames.push_back( words.front() );
	words.pop();
	this->paraFieldTypes.push_back( words.front() );
	words.pop();
    }
DONE:
    return;
}




//
//	readNextLine
//
//	Read the next line of the file, don't do anything with it.
void	MoeReadFile::readNextLine()
{
    char	*cPPut;
    char	c;
    int	iSize;
    if ( this->fMoe->eof() ) {
	this->done = true;
	return;
    }
    cPPut = this->nextLine;
    iSize = sizeof(this->nextLine)-1;
    while ( !this->fMoe->eof() && (c=this->fMoe->get()) && iSize>0 ) 
    {
	if ( c >= ' ' ) {
	    *cPPut = c;
	    cPPut++;
	    iSize--;
	}
	if ( c== 10 ) break;
    }
    *cPPut = '\0';
    if ( this->fMoe->eof() ) {
	this->done = true;
    }
}




//
//	openFileName
//
//	Open the file and read the first line
//
	void	MoeReadFile::openFileName( const char*	sFileName  )
{
    std::ifstream	*fIn;
	struct stat	buf;

	// Create an ifstream and open it

	if ( stat(sFileName,&buf)==-1 ) {
	    string errStr;
	    errStr = sFileName;
	    SIMPLE_ERROR("Could not open MOE file: {}" , errStr );
	}

	fIn = new std::ifstream();
	fIn->open(sFileName);

	this->fileName = sFileName;
	this->fMoe = fIn;
	this->mustDestroy_fMoe = true;
	this->done = false;
	this->readNextLine();
	this->parseHeader();
    }

//
//	openStream
//
//	Create the MoeReadFile from the istream
//
    void	MoeReadFile::openStream( std::istream*	fIn )
    {
	this->fileName = "";
	this->fMoe = fIn;
	this->mustDestroy_fMoe= false;
	this->done = false;
	this->readNextLine();
	this->parseHeader();
    }



//
//	readParaTable
//
//	Fill an empty ParaTable with the data read from the file
    void	MoeReadFile::readParaTable( ParaTable* table )
    {
	QueueString		words;
	char			caStr[MAXWORD];
	ParaInteger*		paraInt;
	ParaArrayInt*		paraArrayInt;
	ParaArrayReal*	paraArrayDbl;
	ParaHexInteger*		paraHex;
	ParaReal*		paraReal;
	ParaToken*		paraToken;
	ParaChar*		paraChar;
	ParaString*		paraString;
	int			row, field;
	string			type;
	string			str;
	char*			cPTemp;
	ParaVector		col;
	int			sz, i;


	// First fill the queue with all the words.
	// this may be inefficient to load all the words
	// for the table.  For now I'll do this.
	do {
	    LOG("readParaTable-- reading next line" );
	    this->readNextLine();
	    LOG("Read line: {}" , this->nextLine  );
	    if ( !this->done && this->nextLine[0] != '#' ) {
		LOG("It's data, splitting it into queue" );
		this->splitLine( words, this->nextLine );
		LOG("Done splitting it" );
	    }
	} while ( !this->done && this->nextLine[0] != '#' );
	LOG("Finished reading data... creating ParaTable" );

	// first Empty the ParaTable
	table->erase();
	for ( field=0; field<this->nextParaTableFieldCount(); field++ ) {
	    LOG("Adding  fieldVector name: {}" , this->nextParaTableFieldName(field).c_str()  );
	    table->addFieldVector( this->nextParaTableFieldName(field), col );
	}

	// Now build the ParaTable

	for ( row=0; row<this->nextParaTableRows(); row++ ) {
	    LOG("Starting row" );
	    for ( field = 0; field < this->nextParaTableFieldCount(); field++ ) {
		LOG("Parsing field: {}" , type.c_str()  );
		type = this->nextParaTableFieldType(field);
		if ( type=="i" ) {
		    str = words.front();
		    words.pop();
		    strncpy( caStr, &(str[0]), str.length() );
		    caStr[str.length()] = '\0';
		    LOG(" i {}" , caStr  );
		    paraInt = new ParaInteger(atoi(caStr));
		    table->getVector(field).push_back(paraInt);
		    continue;
		}
		if ( type.substr(0,2)=="i=" ) {
		    string num;
		    num = type.substr(2,type.length()-2);
		    strncpy( caStr, num.c_str(), num.length() );
		    caStr[num.length()] = '\0';
		    LOG(" i= {}" , caStr  );
		    paraInt = new ParaInteger(atoi(caStr));
		    table->getVector(field).push_back(paraInt);
		    continue;
		}
		if ( type.substr(0,2)=="r=" ) {
		    string num;
		    num = type.substr(2,type.length()-2);
		    strncpy( caStr, num.c_str(), num.length() );
		    caStr[num.length()] = '\0';
		    LOG(" r= {}" , caStr  );
		    paraReal = new ParaReal(atof(caStr));
		    table->getVector(field).push_back(paraReal);
		    continue;
		}
		if ( type=="i*" ) {
		    str = words.front();
		    words.pop();
		    strncpy( caStr, &(str[0]), str.length() );
		    caStr[str.length()] = '\0';
		    paraArrayInt = new ParaArrayInt();
		    sz = atoi(caStr);
		    LOG(" i*= size = {}" , caStr  );
		    int ddd;
		    for ( i=0;i<sz;i++) {
			str = words.front();
			words.pop();
			ddd = atoi(str.c_str());
			LOG(" i*= data[{}] = {}" , i  , str.c_str()  );
			paraArrayInt->push_back(ddd);
		    }
		    table->getVector(field).push_back(paraArrayInt);
		    continue;
		}
		if ( type=="r*" ) {
		    str = words.front();
		    words.pop();
		    strncpy( caStr, &(str[0]), str.length() );
		    caStr[str.length()] = '\0';
		    paraArrayDbl = new ParaArrayReal();
		    sz = atoi(caStr);
		    LOG(" r*= size = {}" , caStr  );
		    double ddd;
		    for ( i=0;i<sz;i++) {
			str = words.front();
			words.pop();
			ddd = atof(str.c_str());
			LOG("   r*= data[{}] = {}" , i  , str.c_str()  );
			paraArrayDbl->push_back(ddd);
		    }
		    table->getVector(field).push_back(paraArrayDbl);
		    continue;
		}
		if ( type=="ix" ) {
		    str = words.front();
		    words.pop();
		    strncpy( caStr, &(str[0]), str.length() );
		    caStr[str.length()] = '\0';
		    LOG("   ix= {}" , caStr  );
		    paraHex = new ParaHexInteger(strtol(caStr,&cPTemp,16));
		    table->getVector(field).push_back(paraHex);
		    continue;
		}
		if ( type=="r" ) {
		    str = words.front();
		    words.pop();
		    strncpy( caStr, &(str[0]), str.length() );
		    caStr[str.length()] = '\0';
		    paraReal = new ParaReal(atof(caStr));
		    LOG("   r = {}" , caStr  );
		    table->getVector(field).push_back(paraReal);
		    continue;
		}
		if ( type=="c" ) {
		    str = words.front();
		    words.pop();
		    paraChar = new ParaChar(str);
		    LOG("   c = {}" , str.c_str()  );
		    table->getVector(field).push_back(paraChar);
		    continue;
		}
		if ( type=="t" ) {
		    str = words.front();
		    words.pop();
		    paraToken = new ParaToken(str);
		    LOG("   t = {}" , str.c_str()  );
		    table->getVector(field).push_back(paraToken);
		    continue;
		}
		if ( type=="s" ) {
		    str = words.front();
		    words.pop();
		    paraString = new ParaString(str);
		    LOG("   s = {}" , str.c_str()  );
		    table->getVector(field).push_back(paraString);
		    continue;
		}
		SIMPLE_ERROR("Unknown MOE file data type({})", type);
	    }
	}
	LOG("Done creating ParaTable" );
	if ( !this->done && this->nextLine[0] == '#' ) {
	    LOG("Parsing header for next ParaTable" );
	    this->parseHeader();
	}
    }





//
//	skipParaTable
//
//	Skip the current paraTable
    void	MoeReadFile::skipParaTable()
    {
	do {
	    this->readNextLine();
	} while ( !this->done && this->nextLine[0] != '#' );
	if ( !this->done && this->nextLine[0] == '#' ) {
	    this->parseHeader();
	}
    }


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



    namespace chem {
	void	moeReadAggregateMoeFile( Aggregate_sp	agg,
					 MoeReadFile	&f,
					 bool		&gotMMTypes );
    };

//
//	moeReadAggregateMoeReadFile
//
//	Read the Aggregate from the MOE file
//
    void	moeReadAggregateMoeFile( Aggregate_sp	agg,
					 MoeReadFile	&f,
					 bool		&gotMMTypes )
    {
	int		i,ir,im,ia,tir,it;
	ParaTable	ptAtoms;
	ParaTable	ptMMTypes;
	ParaTable	ptFixed;
	ParaTable	ptResidues;
	ParaTable	ptMolecules;
	ParaTable	ptBonds;	//!< Read bonds with undefined order
	bool		bGotBonds, bGotBondsO1, bGotBondsO2, bGotBondsO3;
	ParaTable	ptBondsO1;	//!< Read bonds with order o=1
	ParaTable	ptBondsO2;	//!< Read bonds with order o=2
	ParaTable	ptBondsO3;	//!< Read bonds with order o=3
	ParaTable	ptRestraints;
	ParaTable	ptForceRSm1;
	ParaTable	ptForceRSp1;
	ParaTable	ptHintLP;
	ParaTable	ptTemp;
	ParaTable	ptCharges;
        GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Atom_O>>);
        gctools::Vec0<Atom_sp>	atoms;
	int		index_ID;
	int		index_aName;
	int		index_aElement;
	int		index_aGeometry;
	int		index_aPosX;
	int		index_aPosY;
	int		index_aPosZ;
	int		index_a;
	int		index_b;
	int		index_rName;
        int             index_rUID;
	int		index_rAtomCount, index_cResidueCount, numAtoms, numResidues;
	int		index_type;
	int		index_atoms;
	int		index_target;
	int		index_weight;
	string		type;
	double		weight;
	map<string,string>	typeLookup;
	double		x,y,z;
	Vector3		pos;
	int		from, to;
	Atom_sp		rsAtom;
	string		name;
	BondOrder	bo;
	ParaArrayInt*	resAtoms;
	ParaArrayReal*	target;
	bool		bGotRestraints, bGotFixed, bGotCharges;
	bool		bGotForceRSm1, bGotForceRSp1, bGotHintLP;
        core::List_sp	aggRestraints = nil<core::T_O>();
	int		aid;
	int		index_mmTypeId, index_mmTypes;
	string		str;


	bGotBonds = false;
	bGotBondsO1 = false;
	bGotBondsO2 = false;
	bGotBondsO3 = false;
	bGotRestraints = false;
	bGotForceRSm1 = false;
	bGotForceRSp1 = false;
	bGotFixed = false;
	bGotHintLP = false;
	bGotCharges = false;
	gotMMTypes = false;
	while ( f.nextParaTableExists() ) {
	    if ( f.nextParaTableName() == "bond" ) {
		if ( f.nextParaTableFieldCount() == 2 ) {
		    LOG("Reading ptBonds" );
		    f.readParaTable(&ptBonds);
		    bGotBonds = true;
		} else if ( f.nextParaTableFieldCount() == 3 ) {
		    LOG("Checking the order" );
		    string order = f.nextParaTableFieldType(2);
		    if ( order=="i=1" ) {
			f.readParaTable(&ptBondsO1);
			bGotBondsO1 = true;
		    } else if ( order=="i=2" ) {
			f.readParaTable(&ptBondsO2);
			bGotBondsO2 = true;
		    } else if ( order=="i=3" ) {
			f.readParaTable(&ptBondsO3);
			bGotBondsO3 = true;
		    } else {
                      SIMPLE_ERROR("Unknown bond order type: {}", f.nextParaTableFieldType(2));
		    }
		} else {
		    SIMPLE_ERROR("Invalid bond entry in moe file");
		}
	    } else if ( f.nextParaTableName() == "attr" &&
			f.nextParaTableFieldName(1) == "aName" ) {
		f.readParaTable(&ptAtoms);
	    } else if ( f.nextParaTableName() == "attr" &&
			f.nextParaTableFieldName(1) == "aFixed" ) {
		LOG("Starting to read ptFixed" );
		f.readParaTable(&ptFixed);
		LOG("There are {} fixed atoms" , ptFixed.getVector(0).size()  );
		bGotFixed = true;
	    } else if ( f.nextParaTableName() == "attr" &&
			f.nextParaTableFieldName(1) == "aCharge" ) {
		LOG("Starting to read ptCharges" );
		f.readParaTable(&ptCharges);
		LOG("There are {} charged atoms" , ptCharges.getVector(0).size()  );
		bGotCharges = true;
	    } else if ( f.nextParaTableName() == "MMTypes" ) {
		LOG("Reading MMTypes" );
		f.readParaTable(&ptMMTypes);
		gotMMTypes = true;
	    } else if ( f.nextParaTableName() == "attr" &&
			f.nextParaTableFieldName(1) == "rAtomCount" ) {
		f.readParaTable(&ptResidues);
	    } else if ( f.nextParaTableName() == "attr" &&
			f.nextParaTableFieldName(1) == "cResidueCount" ) {
		// PROBLEM reading the next table
		f.readParaTable(&ptMolecules);
	    } else if ( f.nextParaTableName() == "restraint" ) {
		bGotRestraints = true;
		f.readParaTable(&ptRestraints);
	    } else if ( f.nextParaTableName() == "attr" &&
			f.nextParaTableFieldName(1) == "aHintLP" ) {
		LOG("Reading aHintLP" );
		bGotHintLP = true;
		f.readParaTable(&ptHintLP);
	    } else if ( f.nextParaTableName() == "attr" &&
			f.nextParaTableFieldName(1) == "aForceRS" ) {
		// Read the aForceRS paratable
		// First figure out which one it should be
		if ( f.nextParaTableFieldType(1)=="i=-1" ) {
		    LOG("Reading aForceRS i=-1 ParaTable" );
		    bGotForceRSm1 = true;
		    f.readParaTable(&ptForceRSm1);
		} else if ( f.nextParaTableFieldType(1)=="i=1" ) {
		    LOG("Reading aForceRS i=1 ParaTable" );
		    bGotForceRSp1 = true;
		    f.readParaTable(&ptForceRSp1);
		} else {
                    std::cerr << "DAMN IT, THEY CHANGED THE WAY CHIRAL RESTRAINTS ARE STORED! " << std::endl;
                    std::cerr << "should be... aForceRS [i=-1/i=1]" << std::endl;
                    std::cerr << "got... aForceRS " << f.nextParaTableFieldType(1) << std::endl;
		    exit(1);
		}
	    } else {
		LOG("Skipping table with name({})" , f.nextParaTableName().c_str()  );
#ifdef	DEBUG_ON //[
		if ( f.nextParaTableFieldCount() > 0 ) {
		    LOG("Skipping table with FieldType(0): {}" , f.nextParaTableFieldType(0).c_str()  );
		}
		if ( f.nextParaTableFieldCount() > 1 ) {
		    LOG("Skipping table with FieldType(1): {}" , f.nextParaTableFieldType(1).c_str()  );
		}
#endif //]
		f.skipParaTable();
	    }
	}

	LOG("Building the molecule" );

	// Build the molecule
	index_ID = ptAtoms.getFieldIndex("ID");
	index_aName = ptAtoms.getFieldIndex("aName");
	index_aElement = ptAtoms.getFieldIndex("aElement");
	index_aGeometry = ptAtoms.getFieldIndex("aGeometry");
	index_aPosX = ptAtoms.getFieldIndex("aPosX");
	index_aPosY = ptAtoms.getFieldIndex("aPosY");
	index_aPosZ = ptAtoms.getFieldIndex("aPosZ");
	index_rAtomCount = ptResidues.getFieldIndex("rAtomCount");
	index_rName = ptResidues.getFieldIndex("rName");
        index_rUID = ptResidues.getFieldIndex("rUID");
	index_cResidueCount = ptMolecules.getFieldIndex("cResidueCount");


	LOG("Creating atoms" );
	// Create the atoms
	ia = 0;
	tir=0;
	for ( im=0; im<ptMolecules.getVector(0).size(); im++ ) {
	    LOG("Creating molecule #{}" , im  );
	    auto  mol  = gctools::GC<Molecule_O>::allocate_with_default_constructor();
	    LOG(" created" );
	    agg->addMolecule(mol);
            mol->setId(agg->contentSize());
	    LOG(" added to aggregate" );
	    numResidues = dynamic_cast<ParaInteger*>(ptMolecules.getVector(index_cResidueCount)[im])->getValue();
	    LOG("  It has {} residues" , numResidues  );
	    for ( ir=0; ir<numResidues; ir++ ) {
		LOG("     Creating residue: {}" , ir  );
		auto  res  = gctools::GC<Residue_O>::allocate_with_default_constructor();
		//res->setTempInt(ir+1);
		mol->addMatter(res);
		LOG("     Added residue to mol" );
		numAtoms = dynamic_cast<ParaInteger*>(ptResidues.getVector(index_rAtomCount)[tir])->getValue();
		LOG("     Got numAtoms" );
		str = dynamic_cast<ParaToken*>(ptResidues.getVector(index_rName)[tir])->getValue();
		res->setName( chemkw_intern(str) );
                int uid = dynamic_cast<ParaInteger*>(ptResidues.getVector(index_rUID)[tir])->getValue();
                res->setId(uid);
                res->setFileSequenceNumber(uid);
		LOG("     Set residue name: {}" , str.c_str()  );
		LOG("        it has {} atoms" , numAtoms  );
		tir++;
		for ( it=0; it<numAtoms; it++ ) {
                  auto  a  = gctools::GC<Atom_O>::allocate_with_default_constructor();
                    //a->setMoeIndex( dynamic_cast<ParaInteger*>(ptAtoms.getVector(index_ID)[ia])->getValue() );
		    a->setName( chemkw_intern(dynamic_cast<ParaToken*>(ptAtoms .getVector(index_aName)[ia])->getValue() ));
		    a->setElementFromString( dynamic_cast<ParaToken*>(ptAtoms .getVector(index_aElement)[ia])->getValue() );
		    a->setHybridizationFromString( dynamic_cast<ParaToken*>(ptAtoms .getVector(index_aGeometry)[ia])->getValue() );
		    x = dynamic_cast<ParaReal*>(ptAtoms .getVector(index_aPosX)[ia])->getValue();
		    y = dynamic_cast<ParaReal*>(ptAtoms .getVector(index_aPosY)[ia])->getValue();
		    z = dynamic_cast<ParaReal*>(ptAtoms .getVector(index_aPosZ)[ia])->getValue();
		    pos.set( x, y, z );
		    a->setPosition(pos);
		    a->resetFlags();
		    atoms.push_back(a);
		    LOG("Created atom: {}" , a->getName().c_str()  );
		    res->addAtom(a);
                    a->setId(res->contentSize());
		    ia++;
		}
	    }
	}

	//
	// If any atoms are fixed then we need to set that AtomFlag
	//

	LOG("Checking if bGotFixed" );
	if ( bGotFixed) {
	    int i_id;
	    int i_fixed;
	    int fixed;
	    LOG("Yup, got {} fixed atoms" , ptFixed.getVector(0).size()  );
	    i_id    = ptFixed.getFieldIndex("ID");
	    i_fixed = ptFixed.getFieldIndex("aFixed");
	    LOG("Processing Fixed i_id={}, i_fixed={}" , i_id , i_fixed  );
	    for ( i=0; i<ptFixed.getVector(0).size(); i++ ) {
		LOG("Processing element {}" , i  );
		aid = dynamic_cast<ParaInteger*>(ptFixed.getVector(i_id)[i])->getValue()-1;
		LOG("Processing element got aid: {} " , aid  );
		fixed = dynamic_cast<ParaInteger*>(ptFixed.getVector(i_fixed)[i])->getValue();
		LOG("Processing element got fixed: {}" , fixed  );

		if ( fixed != 0 ) {
		    atoms[aid]->modifyFlags(ATOMFLAG_ON,ATOM_FIXED);
		} else {
		    atoms[aid]->modifyFlags(ATOMFLAG_OFF,ATOM_FIXED);
		}
	    }
	    LOG("Done processing aFixed" );
	}

	//
	// If any atoms have charges then we need to set that
	//

	LOG("Checking if bGotCharges" );
	if ( bGotCharges) {
	    int i_id;
	    int i_charge;
	    double charge;
	    LOG("Yup, got {} fixed atoms" , ptCharges.getVector(0).size()  );
	    i_id    = ptCharges.getFieldIndex("ID");
	    i_charge = ptCharges.getFieldIndex("aCharge");
	    LOG("Processing Charges i_id={}, i_charge={}" , i_id , i_charge  );
	    for ( i=0; i<ptCharges.getVector(0).size(); i++ ) {
		LOG("Processing element {}" , i  );
		aid = dynamic_cast<ParaInteger*>(ptCharges.getVector(i_id)[i])->getValue()-1;
		LOG("Processing element got aid: {} " , aid  );
		charge = dynamic_cast<ParaReal*>(ptCharges.getVector(i_charge)[i])->getValue();
		LOG("Processing element got fixed: {}" , fixed  );

		atoms[aid]->setCharge(charge);
	    }
	    LOG("Done processing aCharges" );
	}



	//
	// if we read MMTypes from the moe file then set them in the
	// moeType atom property
	//
	LOG("Checking if gotMMTypes" );
	if ( gotMMTypes ) {
	    LOG("Yup, got MMTypes" );
	    index_mmTypeId = ptMMTypes.getFieldIndex("ID");
	    index_mmTypes = ptMMTypes.getFieldIndex("aMMType");
	    LOG("Processing MMTypes mmTypeId={}, mmTypes={}" , index_mmTypeId , index_mmTypes  );
	    for ( i=0; i<ptMMTypes.getVector(0).size(); i++ ) {
		LOG("Processing element {}" , i  );
		aid = dynamic_cast<ParaInteger*>(ptMMTypes.getVector(index_mmTypeId)[i])->getValue()-1;
		LOG("Processing element got aid: {} " , aid  );
		type = dynamic_cast<ParaToken*>(ptMMTypes.getVector(index_mmTypes)[i])->getValue();
		LOG("Processing element got type: {}" , type.c_str()  );

                IMPLEMENT_MEF("setMoeType");
//		atoms[aid]->setMoeType(chemkw_intern(type));
	    }
	    LOG("Done processing MMTypes" );
	}


	{
	    if ( bGotBonds ) { 
		index_a = ptBonds.getFieldIndex("a");
		index_b = ptBonds.getFieldIndex("b");
		LOG("Got a and b ptBonds.getVector(0).size() = {}" , ptBonds.getVector(0).size()  );
		for ( i=0; i<ptBonds.getVector(0).size(); i++ ) {
		    from = dynamic_cast<ParaInteger*>(ptBonds.getVector(index_a)[i])->getValue()-1;
		    to = dynamic_cast<ParaInteger*>(ptBonds.getVector(index_b)[i])->getValue()-1;
		    LOG("Read bond between atom indexes {}-{}" , from+1 , to+1  );
		    if ( atoms[from]->getHybridization()==hybridization_sp2 &&
			 atoms[to]->getHybridization()==hybridization_sp2 ) 
		    {
			bo = doubleBond;
			atoms[from]->bondTo( atoms[to], bo );
#if 0
			ASSERTP(!atoms[from]->invalid(),
				      "Bond I just formed is invalid atom:"+atoms[from]->getName());
			ASSERTP(!atoms[to]->invalid(),
				      "Bond I just formed is invalid atom:"+atoms[to]->getName());
#endif
		    } else {
			bo = singleBond;
			atoms[from]->bondTo( atoms[to], bo );
#if 0
			ASSERTP(!atoms[from]->invalid(),
				      "Bond I just formed is invalid atom:"+atoms[from]->getName());
			ASSERTP(!atoms[to]->invalid(),
				      "Bond I just formed is invalid atom:"+atoms[to]->getName());
#endif
		    }
		    LOG("Installed bond between atoms: {} - {}" , atoms[from]->description().c_str() , atoms[to]->description().c_str()  );
		}
	    }
	    if ( bGotBondsO1 ) { 
		index_a = ptBondsO1.getFieldIndex("a");
		index_b = ptBondsO1.getFieldIndex("b");
		LOG("Got a and b ptBondsO1.getVector(0).size() = {}" , ptBondsO1.getVector(0).size()  );
		for ( i=0; i<ptBondsO1.getVector(0).size(); i++ ) {
		    from = dynamic_cast<ParaInteger*>(ptBondsO1.getVector(index_a)[i])->getValue()-1;
		    to = dynamic_cast<ParaInteger*>(ptBondsO1.getVector(index_b)[i])->getValue()-1;
		    LOG("Read bond between atom indexes {}-{}" , from+1 , to+1  );
		    bo = singleBond;
		    atoms[from]->bondTo( atoms[to], bo );
		    ASSERT(!atoms[from]->invalid()); // "Bond I just formed is invalid atom:"+atoms[from]->getName());
		    ASSERT(!atoms[to]->invalid()); // "Bond I just formed is invalid atom:"+atoms[to]->getName());
		    LOG("Installed bond between atoms: {} - {}" , atoms[from]->description().c_str() , atoms[to]->description().c_str()  );
		}
	    }
	    if ( bGotBondsO2 ) { 
		index_a = ptBondsO2.getFieldIndex("a");
		index_b = ptBondsO2.getFieldIndex("b");
		LOG("Got a and b ptBondsO2.getVector(0).size() = {}" , ptBondsO2.getVector(0).size()  );
		for ( i=0; i<ptBondsO2.getVector(0).size(); i++ ) {
		    from = dynamic_cast<ParaInteger*>(ptBondsO2.getVector(index_a)[i])->getValue()-1;
		    to = dynamic_cast<ParaInteger*>(ptBondsO2.getVector(index_b)[i])->getValue()-1;
		    LOG("Read bond between atom indexes {}-{}" , from+1 , to+1  );
		    bo = doubleBond;
		    atoms[from]->bondTo( atoms[to], bo );
		    ASSERT(!atoms[from]->invalid()); // "Bond I just formed is invalid atom:"+atoms[from]->getName());
		    ASSERT(!atoms[to]->invalid()); // "Bond I just formed is invalid atom:"+atoms[to]->getName());
		    LOG("Installed bond between atoms: {} - {}" , atoms[from]->description().c_str() , atoms[to]->description().c_str()  );
		}
	    }
	    if ( bGotBondsO3 ) { 
		index_a = ptBondsO3.getFieldIndex("a");
		index_b = ptBondsO3.getFieldIndex("b");
		LOG("Got a and b ptBondsO3.getVector(0).size() = {}" , ptBondsO3.getVector(0).size()  );
		for ( i=0; i<ptBondsO3.getVector(0).size(); i++ ) {
		    from = dynamic_cast<ParaInteger*>(ptBondsO3.getVector(index_a)[i])->getValue()-1;
		    to = dynamic_cast<ParaInteger*>(ptBondsO3.getVector(index_b)[i])->getValue()-1;
		    LOG("Read bond between atom indexes {}-{}" , from+1 , to+1  );
		    bo = tripleBond;
		    atoms[from]->bondTo( atoms[to], bo );
		    ASSERT(!atoms[from]->invalid()); // "Bond I just formed is invalid atom:"+atoms[from]->getName());
		    ASSERT(!atoms[to]->invalid()); // "Bond I just formed is invalid atom:"+atoms[to]->getName());
		    LOG("Installed bond between atoms: {} - {}" , atoms[from]->description().c_str() , atoms[to]->description().c_str()  );
		}
	    }
	}

	LOG("Setting RS force constants" );

	if ( bGotForceRSm1 ) {
	    aggRestraints = agg->allRestraints();
	    IMPLEMENT_ME(); // Handle new restraints
	    index_ID = ptForceRSm1.getFieldIndex("ID");
	    for (i=0; i<ptForceRSm1.getVector(0).size();i++ ) {
		aid = ( dynamic_cast<ParaInteger*>(ptForceRSm1.getVector(index_ID)[i])->getValue()-1 );
		rsAtom = atoms[aid];
		RestraintChiral_sp restrainChiral = RestraintChiral_O::create();
                restrainChiral->setAtomA(rsAtom);
		restrainChiral->setChirality(-1);
		aggRestraints = core::Cons_O::create(restrainChiral,aggRestraints);
	    }
	}
	if ( bGotForceRSp1 ) {
	    aggRestraints = agg->allRestraints();
	    IMPLEMENT_ME(); // handle new restraints
	    index_ID = ptForceRSp1.getFieldIndex("ID");
	    for (i=0; i<ptForceRSp1.getVector(0).size();i++ ) {
		aid = ( dynamic_cast<ParaInteger*>(ptForceRSp1.getVector(index_ID)[i])->getValue()-1 );
		rsAtom = atoms[aid];
		RestraintChiral_sp restrainChiral = RestraintChiral_O::create();
                restrainChiral->setAtomA(rsAtom);
		restrainChiral->setChirality(1);
		aggRestraints = core::Cons_O::create(restrainChiral,aggRestraints);
	    }
	}

	if ( bGotHintLP ) {
	    LOG("Setting HintLP attribute" );
	    index_ID = ptHintLP.getFieldIndex("ID");
	    for ( i=0; i<ptHintLP.getVector(0).size();i++ ) {
		aid = (dynamic_cast<ParaInteger*>(ptHintLP.getVector(index_ID)[i])->getValue()-1 );
		Atom_sp lpAtom = atoms[aid];
		LOG("Setting HintLP=True for atom({})" , lpAtom->getName().c_str() );
		lpAtom->setHintLP(true);
	    }
	}


	if ( bGotRestraints ) {
	    aggRestraints = agg->allRestraints();
	    IMPLEMENT_ME(); // handle new restraints
	    index_type = ptRestraints.getFieldIndex("type");
	    index_atoms = ptRestraints.getFieldIndex("atoms");
	    index_target = ptRestraints.getFieldIndex("target");
	    index_weight = ptRestraints.getFieldIndex("weight");
	    for ( i=0; i<ptRestraints.getVector(0).size(); i++ ) {
		type = dynamic_cast<ParaToken*>
		    (ptRestraints.getVector(index_type)[i])->getValue();
		resAtoms = dynamic_cast<ParaArrayInt*>
		    (ptRestraints.getVector(index_atoms)[i]);
		target = dynamic_cast<ParaArrayReal*>
		    (ptRestraints.getVector(index_target)[i]);
		weight = dynamic_cast<ParaReal*>
		    (ptRestraints.getVector(index_weight)[i])->getValue();

	    	// Process restraint here
		if ( type == "dihedral" ) 
		{
		    RestraintDihedral_sp restraint = RestraintDihedral_O::create();
                    restraint->setAtomA(atoms[(*resAtoms)[0]-1]);
		    restraint->setAtomB(atoms[(*resAtoms)[1]-1]);
		    restraint->setAtomC(atoms[(*resAtoms)[2]-1]);
		    restraint->setAtomD(atoms[(*resAtoms)[3]-1]);
		    restraint->setParameters( (*target)[0], weight );
		    aggRestraints = core::Cons_O::create(restraint,aggRestraints);
		} else if ( type == "angle" ) 
		{
		    RestraintAngle_sp restraint = RestraintAngle_O::create();
		    restraint->setAtomA(atoms[(*resAtoms)[0]-1]);
		    restraint->setAtomB(atoms[(*resAtoms)[1]-1]);
		    restraint->setAtomC(atoms[(*resAtoms)[2]-1]);
		    restraint->setParameters( (*target)[0], (*target)[1], weight );
		    aggRestraints = core::Cons_O::create(restraint,aggRestraints);
		} else if ( type == "distance" ) 
		{
                  SIMPLE_ERROR("Add support for distance restraints");
#if 0
                  RestraintDistance_sp restraint = RestraintDistance_O::create();
		    restraint->setAtomA(atoms[(*resAtoms)[0]-1]);
		    restraint->setAtomB(atoms[(*resAtoms)[1]-1]);
		    restraint->setParameters( (*target)[0], (*target)[1], weight );
		    aggRestraints = core::Cons_O::create(restraint,aggRestraints);
#endif
		} else {
                    std::cerr << "Unknown restraint type: "<< type;
		}
		delete resAtoms;
		delete target;
	    }
	}

	LOG("Done" );
    }


    Aggregate_sp	moeReadAggregateWithAtomTypes(const string& nm)
    {
      auto  agg  = gctools::GC<Aggregate_O>::allocate_with_default_constructor();
	moeReadAggregateWithAtomTypesFromFileName(agg,nm.c_str());
	return agg;
    }


    void	moeReadAggregateWithAtomTypesFromFileName(Aggregate_sp	agg,
							  const char*	sFileName )
    {
	MoeReadFile	f;
	Atom_sp		a;
	Loop		l;
	bool		gotMMTypes;

	LOG("About to open moe file: {}" , sFileName  );
	f.openFileName(sFileName);
	LOG("Opened moe file: {}" , sFileName  );
	moeReadAggregateMoeFile( agg, f, gotMMTypes );

	LOG("Assigning types" );
	if ( gotMMTypes ) {
	    l.loopTopGoal( agg, ATOMS );
	    while ( l.advanceLoopAndProcess() ) {
		a = l.getAtom();
                IMPLEMENT_MEF("setMoeType");
#if 0
		if ( a->getMoeType().notnilp() ) {
		    a->setType(a->getMoeType() );
//		setType( ps, a, a->getMoeType() );
		}
#endif
	    }
	} else {
	    SIMPLE_ERROR("moeReadAggregateWithAtomTypesFromFileName-- there were no types in the file");
	}
	LOG("Done assigning types" );
    }


//
//	moeReadAggregateFromFileName
//
//	Read the Aggregate from the MOE file.
    void	moeReadAggregateFromFileName( Aggregate_sp	agg,
					      const char	*sFileName )
    {
	MoeReadFile	f;
	bool		gotMMTypes;

	LOG("About to open moe file: {}" , sFileName  );
	f.openFileName(sFileName);
	LOG("Opened moe file: {}" , sFileName  );
	moeReadAggregateMoeFile( agg, f, gotMMTypes );
    }


//
//	moeReadAggregate
//
//	Read the Aggregate from the MOE file.
DOCGROUP(cando);
CL_DEFUN  Aggregate_sp chem__moeReadAggregate(const string& name )
{
    MoeReadFile	f;
	bool		gotMMTypes;
	LOG("About to open moe file: {}" , name.c_str() );
	f.openFileName(name.c_str());
	LOG("Opened moe file: {}" , name.c_str() );
	auto  agg  = gctools::GC<Aggregate_O>::allocate_with_default_constructor();
	moeReadAggregateMoeFile( agg, f, gotMMTypes );
	return agg;
    }


//
//	moeReadAggregateFromStream
//
//	Read the Aggregate from the MOE file.
    void	moeReadAggregateFromStream( Aggregate_sp	agg,
					    std::istream		*is)
    {
	MoeReadFile	f;
	bool		gotMMTypes;
	f.openStream(is);
	moeReadAggregateMoeFile( agg, f, gotMMTypes );
    }






//
//	moeWriteAggregateStream
//
//	Write the atoms in the Aggregate to the MOE file
void	moeWriteAggregateStream( Aggregate_sp agg, std::ostream& fout)
    {
	Loop		la,lr,l;
	Atom_sp		a;
	Atom_sp		b;
	Residue_sp	r;
	int		count;
	int		bondCount;
	int		i, resAtomCount;
	BondOrder	order;
        gctools::Vec0<Atom_sp>	atoms;
        gctools::Vec0<Atom_sp>::iterator	ai;
	vector<int>	residueAtomCount;
	vector<int>::iterator	vii;
        gctools::Vec0<MatterName>	residueNames;
	int		resOffset, chainOffset;
	int		fixedCount, hintLPCount;
	Restraint_sp	restraint;
	int		resId;
	// Give each atom a unique ID
	count = 0;
	resId = 1;
	lr.loopTopGoal( agg, RESIDUES );
	{
	    while ( lr.advanceLoopAndProcess() ) {
		r = lr.getResidue();
		
		la.loopTopGoal( r, ATOMS );
		resAtomCount = 0;
		while ( la.advanceLoopAndProcess() ) {
		    a = la.getAtom();
		    count++;
		    resAtomCount++;
		    a->setTempInt(count);
		    LOG("Assigning atom: {} index: {}" , a.get() , count  );
		    atoms.push_back(a);
		};
		residueAtomCount.push_back(resAtomCount);
		residueNames.push_back(r->getPdbName());  // PDBNAME
		//r->setTempInt(resId);
		resId++;
	    }
	}
	chainOffset = count+residueAtomCount.size()+1;
	resOffset = count+1;
	bondCount = 0;
	LOG("Counting bonds" );
	l.loopTopGoal( Matter_sp(agg), BONDS );
	while ( l.advanceLoopAndProcess() ) {
	    LOG("Looking at bond#{}" , bondCount );
	    bondCount++;
	};
	LOG("Writing moe file" );
	fout << "#moe 2001.01" << std::endl;
	fout << "#system 7 tag t value *" << std::endl;
	fout << "BoxEnable i* 3 0 0 0 BoxSize i* 3 10 10 10 BoxShape i* 3 90 90 90" << std::endl;
	fout << "ViewOrientationY r* 3 -0.6040875 0.7739065 0.1901233 ViewOrientationZ r*" << std::endl;
	fout << "3 0.7362871 0.4507356 0.5046968 ViewLookAt r* 3 -0.4651223 0.3664901" << std::endl;
	fout << "-0.2553613 ViewScale r 21.1468887" << std::endl;
	fout << "#molecule 3 tag t value *" << std::endl;
	fout << "atoms i " << count << " residues i "
	     << residueAtomCount.size() << " chains i 1" << std::endl;
	fout << "#bond " << bondCount << " a i b i" << std::endl;
	l.loopTopGoal( Matter_sp(agg), BONDS );
	while ( l.advanceLoopAndProcess() ) {
	    l.getBond( &a, &b, &order );
	    fout << a->getTempInt() << " " << b->getTempInt() << std::endl;
	};
	fout << "#attr " << count << " ID i aName t aElement t aGeometry t aPosX r aPosY r aPosZ r" << std::endl;
	fixedCount = 0;
	hintLPCount = 0;
	LOG("Writing atoms" );
	for ( ai=atoms.begin(); ai!=atoms.end(); ai++ ) {
	    a = (*ai);
	    fout << a->getTempInt() << " ";
	    fout << a->getName() << " ";
	    fout << a->getElement() << " ";
	    // Sometimes hybridization isn't defined.
	    // But it must be present in the output file
	    if ( a->getHybridization() != hybridization_undefined)
	    {
		fout << a->getHybridizationAsString() << " ";
	    } else
	    {
		fout << "sp3 ";
	    }
	    fout << a->getPosition().getX() << " ";
	    fout << a->getPosition().getY() << " ";
	    fout << a->getPosition().getZ() << std::endl;
	    if ( a->getFlags()&ATOM_FIXED ) fixedCount++;
	    if ( a->getHintLP() ) hintLPCount++;
	}

	// Handle Fixed atoms
	if ( fixedCount > 0 ) {
	    fout << "#attr " << fixedCount << " ID i aFixed i=1" << std::endl;
	    for ( ai=atoms.begin(); ai!=atoms.end(); ai++ ) {
		a = (*ai);
		if ( a->getFlags()&ATOM_FIXED ) {
		    fout << a->getTempInt();
		    if ((ai-atoms.begin())%16!=15 ) fout << " ";
		    else fout << std::endl;
		}
	    }
	    fout << std::endl;
	}

	// Handle HintLP atoms
	if ( hintLPCount > 0 ) {
	    fout << "#attr " << hintLPCount << " ID i aHintLP i=1" << std::endl;
	    for ( ai=atoms.begin(); ai!=atoms.end(); ai++ ) {
		a = (*ai);
		if ( a->getHintLP() ) {
		    fout << a->getTempInt();
		    if ((ai-atoms.begin())%16!=15 ) fout << " ";
		    else fout << std::endl;
		}
	    }
	    fout << std::endl;
	}


	fout << "#attr 1 ID i cResidueCount i cName t cHeader t cTag t cColorBy t cRGB ix" << std::endl;
	fout << chainOffset << " " << residueAtomCount.size() << " $ $ $ c:number 0" << std::endl;
	fout << "#attr "<<residueAtomCount.size()<<" ID i rAtomCount i rName t rUID i rINS c rPos i rType t" << std::endl;
	i = 1;
	for ( vii=residueAtomCount.begin();vii!=residueAtomCount.end();vii++ ){
	    fout << resOffset+i-1 << " " << (*vii) << " " ;
	    if ( residueNames[i-1].nilp() ) {
		fout << "$";
	    } else {
		fout << residueNames[i-1];
	    }
	    fout << " 0 ! 0 none" << std::endl;
	    i++;
	}
#if 0
	fout << "#attr 1 ID i rRGB ix rColorBy t rSegment t" << std::endl;
	fout << count+1 << " 0 r:rgb none" << std::endl;
	fout << "#attr 1 ID i rRibbonWidth r rRibbonHeight r rRibbonColorBy t rRibbonMode t" << std::endl;
	fout << count+1 << " 0 0 r:aseg none" << std::endl;
#endif

	// Handle restraints properly
	//
	IMPLEMENT_ME();	

#if 0
	// Process restraints here
	// First count up all of the different types

	restTotal = agg->getRestraints()->size();
	rsm1Count = 0;
	rsp1Count = 0;
	nonRSCount = 0;
	if ( restTotal>0 ) {
	    for (i=0;i<restTotal;i++) {
		restraint = agg->getRestraints()->getRestraintIndex(i);
		if ( restraint->getType()==RESTRAINT_CHIRALITY ) 
		{
		    if ( restraint->getChirality() > 0 ) rsp1Count++;
		    if ( restraint->getChirality() < 0 ) rsm1Count++;
		} else nonRSCount++;
	    }
	}
	// Handle ForceRS+1 restraints
	if ( rsp1Count > 0 ) {
	    fout << "#attr " << rsp1Count << " ID i aForceRS i=1" << std::endl;
	    for (i=0;i<restTotal;i++ ) {
		restraint = agg->getRestraints()->getRestraintIndex(i);
		if ( restraint->getType() == RESTRAINT_CHIRALITY ) {
		    if ( restraint->getChirality() > 0 ) {
			fout << restraint->getAtomA()->getTempInt() << std::endl;
		    }
		}
	    }
	}
	// Handle ForceRS-1 restraints
	if ( rsm1Count > 0 ) {
	    fout << "#attr " << rsm1Count << " ID i aForceRS i=-1" << std::endl;
	    for (i=0;i<restTotal;i++ ) {
		restraint = agg->getRestraints()->getRestraintIndex(i);
		if ( restraint->getType() == RESTRAINT_CHIRALITY ) {
		    if ( restraint->getChirality() < 0 ) {
			fout << restraint->getAtomA()->getTempInt() << std::endl;
		    }
		}
	    }
	}

	fout << "#endmolecule"<<std::endl;

	// Handle non-chiral restraints
	if ( nonRSCount> 0 ) {
	    fout << "#restraint "<<nonRSCount <<" type t atoms i* target r* weight r" << std::endl;
	    for ( i=0; i<restTotal; i++ ) {
		restraint = agg->getRestraints()->getRestraintIndex(i);
		if ( restraint->getType() != RESTRAINT_CHIRALITY ) {
		    _moeWriteRestraint(fout,restraint,agg);
		}
	    }
	}
#endif

	fout << "#endsystem" << std::endl;
    }

//
//	moeWriteAggregateFileName
//
//	Write the atoms in the Aggregate to the MOE file
    void	moeWriteAggregateFileName( Aggregate_sp agg, const char* sFileName )
    {
        std::ofstream	fout;
	fout.open(sFileName);
	moeWriteAggregateStream(agg,fout);
	fout.close();
    }

//
//	moeWriteAggregateReturnAsString
//
//	Return the aggregate in MOE format as a string
//
    string	moeWriteAggregateReturnAsString( Aggregate_sp agg )
    {
	stringstream	fout;
	fout.str("");
	moeWriteAggregateStream(agg,fout);
	return fout.str();
    }



void	moeRead( Aggregate_sp agg, const string& sFileName )
    {
	LOG("Reading from file: {}" , sFileName.c_str()  );
	moeReadAggregateFromFileName(agg, sFileName.c_str() );
	LOG("moeRead done" );
    };


    void	moeWrite( Aggregate_sp a, const string& sFileName ) 
    {
	moeWriteAggregateFileName( a, sFileName.c_str());
    };


};
