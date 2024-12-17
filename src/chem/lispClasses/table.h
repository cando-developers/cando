/*
    File: table.h
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

#ifndef	Table_H //[
#define Table_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
//#include "candoBase/render.fwd.h"

#include "chemPackage.h"


namespace chem {
SMART(Table);

SMART(TableEntry );
class  TableEntry_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,TableEntry_O,"TableEntry");
#if INIT_TO_FACTORIES
 public:
    static TableEntry_sp make();
#else
    DECLARE_INIT();
#endif
public:
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
    gctools::Vec0<core::T_sp>		_Fields;

public:
    static TableEntry_sp createWithFields(uint numberOfFields,Table_sp table, core::LispPtr);

public:
		/*! Return the table (owner) of this entry.
		 */
	Table_sp getTable() const;
	Table_sp getTableNotConst() { return this->getTable();};

		/*! Read the value associated with the field.
		 * Throw exception if it doesn't exist.
		 */
	core::T_sp read(core::Symbol_sp key);

		/*! Put the value associated with the field.
		 * Throw exception if it doesn't exist.
		 */
	void write(core::Symbol_sp key, core::T_sp val);

		/*! Get the value of the field at (fieldIndex)
		 */
	core::T_sp getField(uint fieldIndex);

		/*! Set the value of the field at (fieldIndex)
		 */
	void putField(uint fieldIndex, core::T_sp obj);

		/*! Append an empty field and return its index
		 */
	uint appendField(core::T_sp obj);

		/*! Insert an empty field just before the (position)
		 */
	uint insertField(uint position, core::T_sp obj);

	bool	canRender() { return true; };
#ifdef RENDER
	candoBase::Render_sp rendered(core::Cons_sp options);
#endif

	string __repr__() const;

	TableEntry_O( const TableEntry_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR( TableEntry_O);
};


SMART(Table );
class Table_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,Table_O,"Table");
#if INIT_TO_FACTORIES
 public:
    static Table_sp make();
#else
    DECLARE_INIT();
#endif
public:
public: // virtual functions inherited from Object
	void	initialize();
//	void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
    gctools::Vec0<core::Symbol_sp>		_FieldSymbols;
    gctools::Vec0<core::Class_sp>	_FieldClasses;
    gctools::Vec0<TableEntry_sp>	_Entries;
private: // Don't serialize
    core::HashTable_sp                      _FieldIndices; // m a p<core::Symbol_sp,uint>	_FieldIndices;
public:
    typedef gctools::Vec0<TableEntry_sp>::iterator	entryIterator;
	entryIterator beginEntry() { return this->_Entries.begin(); };
	entryIterator endEntry() { return this->_Entries.end(); };
public:
		/*! Return the index of the field with the fieldSymbol.
		 * If fieldSymbol is nil then return UndefinedUnsignedInt.
		 */
	uint indexOfField(core::Symbol_sp fieldSymbol );
	core::Symbol_sp fieldSymbol(uint idx);
public:
		/*! Return elements as a cons
		 */
	core::Cons_sp asCons();

		/*! Return true if the field exists
		 */
	bool hasField(core::Symbol_sp fieldSymbol);

		/*! Return the Class of the field
		 */
	core::Class_sp fieldClass(core::Symbol_sp fieldSymbol);

		/*! Create a field and append it to the list of fields
		 */
	void appendField(core::Symbol_sp field, core::Class_sp fieldClass );
		/*! Create a field and insert it before (position).
		 * If (position) isNil then insert it as the first field.
		 */
	void insertField(core::Symbol_sp positionField, core::Symbol_sp field, core::Class_sp fieldClass );
		/*! Write the value of the field at the indicated entry
		 */
	void writeField(uint idx, core::Symbol_sp field, core::T_sp value );
		/*! Write the pointed to by symbols to the fields of
		 * the table at the indicated index.
		 * If the Cons contains fieldSymbol/value pairs to write to the
		 * entry. If a field is unknown then an exception is thrown.
		 */
	void write(uint idx, core::Cons_sp values );

		/*! Read the value of the field at the indicated entry
		 */
	core::T_sp readField(uint idx, core::Symbol_sp field );
		/*! Return a binder that maps field symbols to 
		 * values
		 */
	TableEntry_sp read(uint idx);
		/*! Append an entry - return the index.
		 */
	TableEntry_sp appendEntry();
		/*! Insert an empty entry before the given index
		 * Return the index.
		 */
		/*! Append an entry and write the fields to it.
		 */
	void appendWrite(core::Cons_sp elements);
		/*! Insert an empty entry before the given index
		 * Return the index.
		 */
	TableEntry_sp insertEntryBefore(uint beforeIdx);

		/*! Return true if the fields and classes match
		 */
	bool fieldsMatch(Table_sp other);

		/*! Merge the other table into us
		 * If there are different numbers of fields or mismatched fields
		 * then throw an exception.
		 */
	void merge(Table_sp table);

		/*! Return the number of entries
		 */
	uint numberOfEntries();

		/*! Return the contents as a string
		 */
	string contentsAsString();

	string __repr__() const;

	bool	canRender() { return true; };
#ifdef RENDER
	candoBase::Render_sp rendered(core::Cons_sp options);
#endif

public:
	Table_O( const Table_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(Table_O);
};


};
#endif //]
