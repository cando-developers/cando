/*
    File: _msmarts_Parser.cc
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
/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         msmarts_parse
#define yylex           msmarts_lex
#define yyerror         msmarts_error
#define yydebug         msmarts_debug
#define yynerrs         msmarts_nerrs


/* Copy the first part of user declarations.  */
#line 69 "msmarts_Parser.yy" /* yacc.c:339  */

#define DEBUG_LEVEL_NONE

#include<iostream>
#include<string>
#include<vector>
#include <istream>

#include <clasp/core/common.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>

#include <cando/chem/chemInfo.h>
#include <cando/chem/hold.h>


using namespace std;
using namespace chem;

// This defines the stream that we are reading from

#include <cando/chem/msmarts_TypeParser.h>
#include <cando/chem/msmarts_ParserParam.h>
#include <cando/chem/msmarts_LocationType.h>

#define	LEXDEBUG	1

#ifdef	LEXDEBUG
#define	LEXPRINT(x) {lisp_LOG( "Token: {}" %(x));}
#define	LEXDPRINT(x) {lisp_LOG(x)};
#else
#define	LEXPRINT(x)
#define	LEXDPRINT(x)
#endif
#define	MAX_CVAL	16

#define MP() (data)

void msmarts_error(YYLTYPE* yyloc, msmarts_SParserParam* data, const char* message );
int	msmarts_lex(YYSTYPE* yylval, YYLTYPE* yylloc, msmarts_SParserParam* data);


#line 117 "_msmarts_Parser.cc" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "_msmarts_Parser.hh".  */
#ifndef YY_MSMARTS_MSMARTS_PARSER_HH_INCLUDED
# define YY_MSMARTS_MSMARTS_PARSER_HH_INCLUDED
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

#endif /* !YY_MSMARTS_MSMARTS_PARSER_HH_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 227 "_msmarts_Parser.cc" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
             && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  49
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   121

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  51
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  12
/* YYNRULES -- Number of rules.  */
#define YYNRULES  74
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  93

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   305

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   206,   206,   211,   212,   216,   223,   227,   231,   235,
     241,   245,   266,   270,   275,   282,   288,   295,   301,   305,
     310,   336,   347,   353,   365,   376,   388,   393,   404,   409,
     421,   433,   444,   456,   470,   482,   486,   497,   501,   510,
     514,   524,   528,   534,   538,   542,   546,   559,   563,   573,
     582,   598,   602,   611,   620,   624,   628,   632,   641,   645,
     649,   658,   662,   666,   672,   676,   680,   689,   766,   776,
     786,   796,   807,   814,   817
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "APEndOfLine", "APBond",
  "APAtomicNumber", "APChiralityAntiClockwise", "APChiralityClockwise",
  "APCloseBracket", "APCloseParenthesis", "APConnectivity", "APDegree",
  "APRingTest", "APResidueTest", "APDollarSign", "APDoubleBond",
  "APAliphatic", "APAromatic", "APOrganicElement",
  "APAromaticOrganicElement", "APInorganicElement",
  "APElectronegativeElement", "APError", "APGroupNumber",
  "APHeavyAtomTotalBond", "APImplicitHCount", "AM1_BCC_x", "AM1_BCC_y",
  "APLonePair", "APNegativeCharge", "APNegativeFormalCharge", "APNumber",
  "APLambda", "APOpenBracket", "APOpenParenthesis", "APOperatorAndHigh",
  "APOperatorAndLow", "APOperatorNot", "APOperatorOr", "APPiBondOrbital",
  "APAromaticPiElectron", "APPositiveCharge", "APPositiveFormalCharge",
  "APRingMembershipCount", "APRingSize", "APTotalHCount",
  "APTotalBondNumber", "APTotalExplicitBondNumber", "APTransitionMetal",
  "APValence", "APWildCard", "$accept", "input", "chemInfo", "chain",
  "bondAtomTest", "acyclicAtomTest", "atomTest", "branch",
  "atomPrimativeTest", "logOp", "recursiveChemInfo", "intNumber", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305
};
# endif

#define YYPACT_NINF -50

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-50)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
      46,   -26,   -50,   -50,   -50,   -50,    71,    11,   -50,   -14,
      23,   -50,    46,    -9,   -50,   -50,    -9,    -9,    -5,    -3,
     -50,   -50,   -50,   -50,   -50,   -50,    -9,    -9,    -9,   -50,
     -50,    -9,   -13,   -50,   -50,    71,   -50,    -9,   -27,    -9,
      -9,    -9,    -9,    -9,   -50,   -50,    15,   -50,   -50,   -50,
     -50,    46,    23,   -50,    23,   -50,    23,    29,   -50,   -50,
     -50,   -50,   -50,   -50,   -50,   -50,   -50,    14,   -50,   -50,
     -50,     3,   -50,   -50,   -50,   -50,   -50,   -50,   -50,    71,
      71,    71,   -50,    36,   -50,   -50,   -50,   -50,   -50,   -50,
     -23,    19,   -50
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       3,     0,    15,    16,    13,    14,     0,     0,     2,    19,
       4,    17,     3,     0,    59,    60,    41,    22,     0,     0,
      31,    32,    62,    63,    64,    44,     0,    45,    28,    24,
      25,    65,    48,    74,    33,     0,    57,     0,    52,    35,
      37,    26,     0,    39,    21,    68,     0,    67,    61,     1,
      18,     0,     0,     5,     8,    11,     9,     0,    58,    42,
      23,    30,    34,    43,    46,    29,    66,    49,    47,    69,
      56,    53,    51,    36,    38,    27,    55,    40,    12,     0,
       0,     0,    10,     0,     6,     7,    73,    50,    54,    70,
      72,    71,    20
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -50,   -50,    37,   -49,   -50,   -50,     1,   -50,   -50,   -33,
      -6,    -7
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     7,     8,    53,    54,     9,    55,    56,    45,    46,
      11,    48
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      47,    10,    69,    83,    33,    84,    58,    85,    12,    59,
      60,    49,    79,    10,    71,    81,    67,    50,    33,    63,
      64,    65,    33,    78,    66,    68,    61,    51,    62,    47,
      70,    72,    73,    74,    75,    76,    77,     1,    86,     2,
       3,     4,     5,    87,    88,    92,    89,    90,    91,    57,
      79,    80,    82,    81,    79,     0,     6,    52,     0,     0,
       1,     0,     2,     3,     4,     5,     0,     0,     0,     0,
       0,     0,     0,    47,    47,    47,    13,    14,    15,     6,
       0,    16,    17,    18,    19,     1,     0,    20,    21,    22,
      23,    24,    25,     0,    26,    27,    28,    29,    30,    31,
      32,     0,    33,    34,     0,     0,     0,     0,    35,     0,
      36,    37,    38,     0,    39,    40,    41,     0,     0,    42,
      43,    44
};

static const yytype_int8 yycheck[] =
{
       6,     0,    35,    52,    31,    54,    13,    56,    34,    16,
      17,     0,    35,    12,    41,    38,    29,    31,    31,    26,
      27,    28,    31,     8,    31,    32,    31,     4,    31,    35,
      37,    38,    39,    40,    41,    42,    43,    14,     9,    16,
      17,    18,    19,    29,    41,     9,    79,    80,    81,    12,
      35,    36,    51,    38,    35,    -1,    33,    34,    -1,    -1,
      14,    -1,    16,    17,    18,    19,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    79,    80,    81,     5,     6,     7,    33,
      -1,    10,    11,    12,    13,    14,    -1,    16,    17,    18,
      19,    20,    21,    -1,    23,    24,    25,    26,    27,    28,
      29,    -1,    31,    32,    -1,    -1,    -1,    -1,    37,    -1,
      39,    40,    41,    -1,    43,    44,    45,    -1,    -1,    48,
      49,    50
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    14,    16,    17,    18,    19,    33,    52,    53,    56,
      57,    61,    34,     5,     6,     7,    10,    11,    12,    13,
      16,    17,    18,    19,    20,    21,    23,    24,    25,    26,
      27,    28,    29,    31,    32,    37,    39,    40,    41,    43,
      44,    45,    48,    49,    50,    59,    60,    61,    62,     0,
      31,     4,    34,    54,    55,    57,    58,    53,    62,    62,
      62,    31,    31,    62,    62,    62,    62,    29,    62,    60,
      62,    41,    62,    62,    62,    62,    62,    62,     8,    35,
      36,    38,    57,    54,    54,    54,     9,    29,    41,    60,
      60,    60,     9
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    51,    52,    53,    53,    53,    54,    54,    54,    54,
      55,    55,    56,    56,    56,    56,    56,    56,    57,    57,
      58,    59,    59,    59,    59,    59,    59,    59,    59,    59,
      59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
      59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
      59,    59,    59,    59,    59,    59,    59,    59,    59,    59,
      59,    59,    59,    59,    59,    59,    59,    59,    60,    60,
      60,    60,    60,    61,    62
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     0,     1,     2,     2,     2,     1,     1,
       2,     1,     3,     1,     1,     1,     1,     1,     2,     1,
       3,     1,     1,     2,     1,     1,     1,     2,     1,     2,
       2,     1,     1,     1,     2,     1,     2,     1,     2,     1,
       2,     1,     2,     2,     1,     1,     2,     2,     1,     2,
       3,     2,     1,     2,     3,     2,     2,     1,     2,     1,
       1,     1,     1,     1,     1,     1,     2,     1,     1,     2,
       3,     3,     3,     4,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (&yylloc, data, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (0)
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL

/* Print *YYLOCP on YYO.  Private, do not rely on its existence. */

YY_ATTRIBUTE_UNUSED
static unsigned
yy_location_print_ (FILE *yyo, YYLTYPE const * const yylocp)
{
  unsigned res = 0;
  int end_col = 0 != yylocp->last_column ? yylocp->last_column - 1 : 0;
  if (0 <= yylocp->first_line)
    {
      res += YYFPRINTF (yyo, "%d", yylocp->first_line);
      if (0 <= yylocp->first_column)
        res += YYFPRINTF (yyo, ".%d", yylocp->first_column);
    }
  if (0 <= yylocp->last_line)
    {
      if (yylocp->first_line < yylocp->last_line)
        {
          res += YYFPRINTF (yyo, "-%d", yylocp->last_line);
          if (0 <= end_col)
            res += YYFPRINTF (yyo, ".%d", end_col);
        }
      else if (0 <= end_col && yylocp->first_column < end_col)
        res += YYFPRINTF (yyo, "-%d", end_col);
    }
  return res;
 }

#  define YY_LOCATION_PRINT(File, Loc)          \
  yy_location_print_ (File, &(Loc))

# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, Location, data); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, msmarts_SParserParam* data)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (yylocationp);
  YYUSE (data);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, msmarts_SParserParam* data)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  YY_LOCATION_PRINT (yyoutput, *yylocationp);
  YYFPRINTF (yyoutput, ": ");
  yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp, data);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, msmarts_SParserParam* data)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                       , &(yylsp[(yyi + 1) - (yynrhs)])                       , data);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, yylsp, Rule, data); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, msmarts_SParserParam* data)
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);
  YYUSE (data);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yytype)
    {
          case 53: /* chemInfo  */
#line 194 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).smartsRoot)!=NULL) delete (((*yyvaluep).smartsRoot));}
#line 1211 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;

    case 54: /* chain  */
#line 196 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).bondListMatchNode)!=NULL) delete (((*yyvaluep).bondListMatchNode)); }
#line 1217 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;

    case 55: /* bondAtomTest  */
#line 200 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).bondTest)!=NULL) delete (((*yyvaluep).bondTest)); }
#line 1223 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;

    case 56: /* acyclicAtomTest  */
#line 201 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).atomOrBondMatchNode)!=NULL) delete (((*yyvaluep).atomOrBondMatchNode)); }
#line 1229 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;

    case 57: /* atomTest  */
#line 198 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).atomOrBondMatchNode)!=NULL) delete (((*yyvaluep).atomOrBondMatchNode)); }
#line 1235 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;

    case 58: /* branch  */
#line 197 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).bondListMatchNode)!=NULL) delete (((*yyvaluep).bondListMatchNode)); }
#line 1241 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;

    case 59: /* atomPrimativeTest  */
#line 199 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).atomOrBondMatchNode)!=NULL) delete (((*yyvaluep).atomOrBondMatchNode)); }
#line 1247 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;

    case 60: /* logOp  */
#line 202 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).logical)!=NULL) delete (((*yyvaluep).logical)); }
#line 1253 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;

    case 61: /* recursiveChemInfo  */
#line 195 "msmarts_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).root)!=NULL) delete (((*yyvaluep).root)); }
#line 1259 "_msmarts_Parser.cc" /* yacc.c:1257  */
        break;


      default:
        break;
    }
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (msmarts_SParserParam* data)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

/* Location data for the lookahead symbol.  */
static YYLTYPE yyloc_default
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  = { 1, 1, 1, 1 }
# endif
;
YYLTYPE yylloc = yyloc_default;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.
       'yyls': related to locations.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls;
    YYLTYPE *yylsp;

    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[3];

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yylsp = yyls = yylsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  yylsp[0] = yylloc;
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;
        YYLTYPE *yyls1 = yyls;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yyls1, yysize * sizeof (*yylsp),
                    &yystacksize);

        yyls = yyls1;
        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
        YYSTACK_RELOCATE (yyls_alloc, yyls);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex (&yylval, &yylloc, data);
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  *++yylsp = yylloc;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location.  */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 207 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("input: chemInfo");
	    data->expression = (yyvsp[0].smartsRoot)->_obj; 
	}
#line 1555 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 3:
#line 211 "msmarts_Parser.yy" /* yacc.c:1661  */
    { (yyval.smartsRoot) = new Hold<chem::SmartsRoot_O>(); }
#line 1561 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 4:
#line 213 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("chemInfo: atomTest");
		    (yyval.smartsRoot) = new Hold<chem::SmartsRoot_O>(chem::SmartsRoot_O::create((yyvsp[0].atomOrBondMatchNode)->_obj)); 
		}
#line 1569 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 5:
#line 217 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("chemInfo:atomTest chain");
		    (yyval.smartsRoot) = new Hold<chem::SmartsRoot_O>(chem::SmartsRoot_O::create((yyvsp[-1].atomOrBondMatchNode)->_obj,(yyvsp[0].bondListMatchNode)->_obj)); 
		}
#line 1577 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 6:
#line 224 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("chain:bondAtomTest chain");
	(yyval.bondListMatchNode) = new Hold<chem::BondListMatchNode_O>(chem::Chain_O::create( (yyvsp[-1].bondTest)->_obj, (yyvsp[0].bondListMatchNode)->_obj)); 
    }
#line 1585 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 7:
#line 228 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("chain: branch chain");
	(yyval.bondListMatchNode) = new Hold<chem::BondListMatchNode_O>(chem::Branch_O::create( (yyvsp[-1].bondListMatchNode)->_obj, (yyvsp[0].bondListMatchNode)->_obj)); 
    }
#line 1593 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 8:
#line 232 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("chain: bondAtomTest");
	(yyval.bondListMatchNode) = new Hold<chem::BondListMatchNode_O>(chem::Chain_O::create( (yyvsp[0].bondTest)->_obj)); 
    }
#line 1601 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 9:
#line 236 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("chain: branch");
	(yyval.bondListMatchNode) = new Hold<chem::BondListMatchNode_O>(chem::Branch_O::create( (yyvsp[0].bondListMatchNode)->_obj)); 
    }
#line 1609 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 10:
#line 242 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("bondAtomTest:APBond atomTest");
			    (yyval.bondTest) = new Hold<chem::BondTest_O>(chem::BondTest_O::create( (yyvsp[-1].benum), (yyvsp[0].atomOrBondMatchNode)->_obj)); 
			}
#line 1617 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 11:
#line 246 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("bondAtomTest: atomTest");
			    (yyval.bondTest) = new Hold<chem::BondTest_O>(chem::BondTest_O::create( chem::SABSingleOrAromaticBond, (yyvsp[0].atomOrBondMatchNode)->_obj)); 
			}
#line 1625 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 12:
#line 267 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("acyclicAtomTest: APOpenBracket logOp APCloseBracket");
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>((yyvsp[-1].logical)->_obj);
	}
#line 1633 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 13:
#line 271 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("acyclicAtomTest: APOrganicElement");
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>((chem::AtomTest_O::create(chem::SAPElement, _lisp->internKeyword((yyvsp[0].eval))))); 
	}
#line 1641 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 14:
#line 276 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("acyclicAtomTest: APAromaticOrganicElement");
	    string elementName = (yyvsp[0].eval);
	    elementName[0] = toupper(elementName[0]);
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>((chem::AtomTest_O::create(chem::SAPAromaticElement, _lisp->internKeyword(elementName)))); 
	}
#line 1651 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 15:
#line 283 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
	    (yyval.atomOrBondMatchNode)=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPAliphatic));
	}
#line 1660 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 16:
#line 289 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
	    (yyval.atomOrBondMatchNode)=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPAromatic));
	}
#line 1669 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 17:
#line 296 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("acyclicAtomTest: recursiveChemInfo");
	    		    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>((yyvsp[0].root)->_obj); 
			}
#line 1677 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 18:
#line 302 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACEF(BF("atomTest: acyclicAtomTest APNumber   // atomTag==(%s)  ") % ((yyvsp[0].cval)));
	(yyval.atomOrBondMatchNode)=new Hold<chem::AtomOrBondMatchNode_O>(chem::TagSet_O::create((yyvsp[-1].atomOrBondMatchNode)->_obj,(yyvsp[0].cval)));
    }
#line 1685 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 19:
#line 306 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("atomTest: acyclicAtomTest");
	(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>((yyvsp[0].atomOrBondMatchNode)->_obj); 
    }
#line 1693 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 20:
#line 310 "msmarts_Parser.yy" /* yacc.c:1661  */
    { (yyval.bondListMatchNode) = new Hold<chem::BondListMatchNode_O>((yyvsp[-1].bondListMatchNode)->_obj); }
#line 1699 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 21:
#line 337 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPWildCard)); 
	}
#line 1708 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 22:
#line 348 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( 
				chem::AtomTest_O::create( chem::SAPDegree, 1)); 
	}
#line 1717 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 23:
#line 354 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( 
			chem::AtomTest_O::create( chem::SAPDegree,(yyvsp[0].ival))); 
	}
#line 1726 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 24:
#line 366 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create( chem::SAPAM1_BCC_x, 0));
	}
#line 1735 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 25:
#line 377 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create( chem::SAPAM1_BCC_y, 0));
	}
#line 1744 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 26:
#line 389 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create( chem::SAPTotalHCount, 1));
	}
#line 1753 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 27:
#line 394 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create( chem::SAPTotalHCount, (yyvsp[0].ival)));
	}
#line 1762 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 28:
#line 405 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( 
			chem::AtomTest_O::create( chem::SAPImplicitHCount, 1)); 
	}
#line 1771 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 29:
#line 410 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
		(yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>(
					chem::AtomTest_O::create( 
						chem::SAPImplicitHCount, (yyvsp[0].ival)));
	}
#line 1781 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 30:
#line 422 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
		(yyval.atomOrBondMatchNode)=new Hold<chem::AtomOrBondMatchNode_O>(
					chem::AtomTest_O::create(chem::SAPRingTest,(yyvsp[0].cval)));
	}
#line 1790 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 31:
#line 434 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
	    (yyval.atomOrBondMatchNode)=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPAliphatic));
	}
#line 1799 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 32:
#line 445 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
	    (yyval.atomOrBondMatchNode)=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPAromatic));
	}
#line 1808 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 33:
#line 457 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
	    string symStr((yyvsp[0].cval));
	    core::Symbol_sp sym = _lisp->intern(symStr);
	    (yyval.atomOrBondMatchNode)=new Hold<chem::AtomOrBondMatchNode_O>(
				chem::AtomTest_O::create(chem::SAPLambda,sym));
	}
#line 1819 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 34:
#line 471 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
		(yyval.atomOrBondMatchNode)=new Hold<chem::AtomOrBondMatchNode_O>(
					chem::AtomTest_O::create(chem::SAPResidueTest,(yyvsp[0].cval)));
	}
#line 1828 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 35:
#line 483 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPRingMembershipCount, 1)); 
	}
#line 1836 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 36:
#line 487 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPRingMembershipCount, (yyvsp[0].ival)));
	}
#line 1844 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 37:
#line 498 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPRingSize, 1)); 
	}
#line 1852 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 38:
#line 502 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPRingSize, (yyvsp[0].ival)));
	}
#line 1860 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 39:
#line 511 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPValence, 1)); 
	}
#line 1868 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 40:
#line 515 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPValence, (yyvsp[0].ival)));
	}
#line 1876 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 41:
#line 525 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPConnectivity, 1)); 
	}
#line 1884 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 42:
#line 529 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPConnectivity, (yyvsp[0].ival)));
	}
#line 1892 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 43:
#line 535 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPGroupNumber, (yyvsp[0].ival))); 
	}
#line 1900 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 44:
#line 539 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPElectronegativeElement)); 
	}
#line 1908 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 45:
#line 543 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPHeavyAtomTotalBond, 1)); 
	}
#line 1916 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 46:
#line 547 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPHeavyAtomTotalBond, (yyvsp[0].ival))); 
	}
#line 1924 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 47:
#line 560 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeCharge, (yyvsp[0].ival)));
	}
#line 1932 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 48:
#line 564 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeCharge, -1));
	}
#line 1940 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 49:
#line 574 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeCharge, -2));
	}
#line 1948 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 50:
#line 583 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPNegativeCharge, -3));
	}
#line 1956 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 51:
#line 599 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveCharge, (yyvsp[0].ival)));
	}
#line 1964 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 52:
#line 603 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveCharge, 1));
	}
#line 1972 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 53:
#line 612 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveCharge, 2));
	}
#line 1980 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 54:
#line 621 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPositiveCharge, 3));
	}
#line 1988 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 55:
#line 625 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPTransitionMetal, (yyvsp[0].ival))); 
	}
#line 1996 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 56:
#line 629 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPAromaticPiElectron , (yyvsp[0].ival)));
	}
#line 2004 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 57:
#line 633 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPPiBondOrbital)); 
	}
#line 2012 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 58:
#line 642 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPAtomicNumber, (yyvsp[0].ival)));
	}
#line 2020 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 59:
#line 646 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPChiralityAntiClockwise));
	}
#line 2028 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 60:
#line 650 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPChiralityClockwise));
	}
#line 2036 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 61:
#line 659 "msmarts_Parser.yy" /* yacc.c:1661  */
    {
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPAtomicMass, (yyvsp[0].ival)));
	}
#line 2044 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 62:
#line 663 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPElement, _lisp->internKeyword((yyvsp[0].eval)) )); 
	}
#line 2052 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 63:
#line 667 "msmarts_Parser.yy" /* yacc.c:1661  */
    { _lisp_BLOCK_TRACE("atomPrimitiveTest: APAromaticOrganicElement");
	    string elementName = (yyvsp[0].eval);
	    elementName[0] = toupper(elementName[0]);
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPAromaticElement, _lisp->internKeyword(elementName) )); 
	}
#line 2062 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 64:
#line 673 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPElement, _lisp->internKeyword((yyvsp[0].eval)) )); 
	}
#line 2070 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 65:
#line 677 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPLonePair, 1)); 
	}
#line 2078 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 66:
#line 681 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>( chem::AtomTest_O::create( chem::SAPLonePair, (yyvsp[0].ival)));
	}
#line 2086 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 67:
#line 690 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.atomOrBondMatchNode) = new Hold<chem::AtomOrBondMatchNode_O>((yyvsp[0].root)->_obj); 
	}
#line 2094 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 68:
#line 767 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.logical) = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logIdentity, (yyvsp[0].atomOrBondMatchNode)->_obj )); 
	}
#line 2103 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 69:
#line 777 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.logical) = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logNot, (yyvsp[0].logical)->_obj )); 
	}
#line 2112 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 70:
#line 787 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.logical) = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logHighPrecedenceAnd, (yyvsp[-2].logical)->_obj, (yyvsp[0].logical)->_obj )); 
	}
#line 2121 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 71:
#line 797 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.logical) = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logOr, (yyvsp[-2].logical)->_obj, (yyvsp[0].logical)->_obj )); 
	}
#line 2130 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 72:
#line 808 "msmarts_Parser.yy" /* yacc.c:1661  */
    { 
	    (yyval.logical) = new Hold<chem::Logical_O>(
	    		chem::Logical_O::create( chem::logLowPrecedenceAnd, (yyvsp[-2].logical)->_obj, (yyvsp[0].logical)->_obj )); 
	}
#line 2139 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 73:
#line 814 "msmarts_Parser.yy" /* yacc.c:1661  */
    { (yyval.root) = new Hold<chem::Root_O>((yyvsp[-1].smartsRoot)->_obj); }
#line 2145 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;

  case 74:
#line 817 "msmarts_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = atoi((yyvsp[0].cval)); }
#line 2151 "_msmarts_Parser.cc" /* yacc.c:1661  */
    break;


#line 2155 "_msmarts_Parser.cc" /* yacc.c:1661  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (&yylloc, data, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (&yylloc, data, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }

  yyerror_range[1] = yylloc;

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, &yylloc, data);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  yyerror_range[1] = yylsp[1-yylen];
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, yylsp, data);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  yyerror_range[2] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the lookahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
  *++yylsp = yyloc;

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (&yylloc, data, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, &yylloc, data);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, yylsp, data);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 822 "msmarts_Parser.yy" /* yacc.c:1906  */


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
	            SIMPLE_ERROR("Exceeded max number of digits in Lambda name");
		}
		yylval->cval[charCount] = ichar;
		charCount++;
	    }
	    ichar = data->lexStream->get();
	    if ( ichar != '>' )
	    {
	        SIMPLE_ERROR("Illegal lambda name must match <[A-Za-z0-9]+> - instead saw[{}]" , ichar );
	    }	
	    yylval->cval[charCount] = '\0';
	    lisp_LOG("APLambda<{}>" , (yylval->cval) );
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
		        SIMPLE_ERROR("Exceeded max number of digits in Number");
		    }
		    yylval->cval[digitCount] = ichar;
		    digitCount++;
	    }
	    yylval->cval[digitCount] = '\0';
	    lisp_LOG("APNumber = ({})" , (yylval->cval) );
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


chem::SmartsRoot_sp smarts_compile(const string& input, stringstream& errorStream)
{
    msmarts_SParserParam p;
    stringstream sin(input);
    p.lisp = lisp;
    p.expression = _Nil<SmartsRoot_O>();
    p.lexStream = &sin;
    p.msmartsErrorStream = &errorStream;
    lisp_LOG("Entering msmarts_parse");
    if ( msmarts_parse(&p) )
    {
	SIMPLE_ERROR("{}" , (p.msmartsErrorStream->str()));
    }
    return p.expression;
}


void msmarts_error(YYLTYPE* yyloc, msmarts_SParserParam* data, const char* message )
{
    *(data->msmartsErrorStream) << "Error: " << message << endl;
}
