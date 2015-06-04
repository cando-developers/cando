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

#ifndef YY_MSMARTS_MSMARTS_PARSER_HPP_INCLUDED
# define YY_MSMARTS_MSMARTS_PARSER_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int msmarts_debug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    APEndOfLine = 258,
    APBond = 259,
    APAtomicNumber = 260,
    APChiralityAntiClockwise = 261,
    APChiralityClockwise = 262,
    APCloseBracket = 263,
    APCloseParenthesis = 264,
    APConnectivity = 265,
    APDegree = 266,
    APRingTest = 267,
    APResidueTest = 268,
    APDollarSign = 269,
    APDoubleBond = 270,
    APAliphatic = 271,
    APAromatic = 272,
    APOrganicElement = 273,
    APAromaticOrganicElement = 274,
    APInorganicElement = 275,
    APElectronegativeElement = 276,
    APError = 277,
    APGroupNumber = 278,
    APHeavyAtomTotalBond = 279,
    APImplicitHCount = 280,
    AM1_BCC_x = 281,
    AM1_BCC_y = 282,
    APLonePair = 283,
    APNegativeCharge = 284,
    APNegativeFormalCharge = 285,
    APNumber = 286,
    APLambda = 287,
    APOpenBracket = 288,
    APOpenParenthesis = 289,
    APOperatorAndHigh = 290,
    APOperatorAndLow = 291,
    APOperatorNot = 292,
    APOperatorOr = 293,
    APPiBondOrbital = 294,
    APAromaticPiElectron = 295,
    APPositiveCharge = 296,
    APPositiveFormalCharge = 297,
    APRingMembershipCount = 298,
    APRingSize = 299,
    APTotalHCount = 300,
    APTotalBondNumber = 301,
    APTotalExplicitBondNumber = 302,
    APTransitionMetal = 303,
    APValence = 304,
    APWildCard = 305
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



int msmarts_parse (msmarts_SParserParam* data);

#endif /* !YY_MSMARTS_MSMARTS_PARSER_HPP_INCLUDED  */
