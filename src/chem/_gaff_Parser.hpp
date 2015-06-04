/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_GAFF_GAFF_PARSER_HPP_INCLUDED
# define YY_GAFF_GAFF_PARSER_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int gaff_debug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    APEndOfLine = 258,
    APElement = 259,
    APElementWildCard = 260,
    APWildCard = 261,
    APOpenParenthesis = 262,
    APCloseParenthesis = 263,
    APOpenBracket = 264,
    APOpenAngleBracket = 265,
    APCloseBracket = 266,
    APCloseAngleBracket = 267,
    APOperatorAnd = 268,
    APOperatorOr = 269,
    APColon = 270,
    APAntechamberARLevel = 271,
    APAntechamberRingMembership = 272,
    APAntechamberNoRingMembership = 273,
    APBond = 274,
    APNumber = 275,
    APTag = 276,
    APError = 277,
    APNotBondedToPrevious = 278,
    APBondedToPrevious = 279,
    APAtomTypeDefinition = 280,
    APOperatorNot = 281
  };
#endif

/* Value type.  */

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif



int gaff_parse (gaff_SParserParam* data);

#endif /* !YY_GAFF_GAFF_PARSER_HPP_INCLUDED  */
