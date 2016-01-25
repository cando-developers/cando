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
#define yyparse         gaff_parse
#define yylex           gaff_lex
#define yyerror         gaff_error
#define yydebug         gaff_debug
#define yynerrs         gaff_nerrs


/* Copy the first part of user declarations.  */
#line 9 "gaff_Parser.yy" /* yacc.c:339  */

#define DEBUG_LEVEL_NONE
#include <string.h>
#include <iostream>
#include <string>
#include <vector>
#include <istream>

#include "core/foundation.h"
#include "aggregate.h"
#include "molecule.h"
#include "residue.h"
#include "atom.h"
#include "hold.h"

#include "chemInfo.h"
#include "gaff_TypeParser.h"
#include "gaff_ParserParam.h"
#include "gaff_LocationType.h"


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




#line 120 "_gaff_Parser.cpp" /* yacc.c:339  */

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
   by #include "_gaff_Parser.hpp".  */
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

/* Copy the second part of user declarations.  */

#line 206 "_gaff_Parser.cpp" /* yacc.c:358  */

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
#define YYFINAL  8
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   73

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  27
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  22
/* YYNRULES -- Number of rules.  */
#define YYNRULES  56
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  88

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   281

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
      25,    26
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   141,   141,   146,   162,   177,   191,   206,   218,   229,
     239,   248,   259,   260,   261,   262,   266,   267,   274,   279,
     280,   284,   285,   289,   290,   294,   295,   299,   300,   307,
     308,   314,   318,   319,   320,   323,   324,   325,   326,   327,
     328,   329,   330,   331,   335,   338,   339,   343,   344,   349,
     355,   363,   364,   367,   368,   372,   373
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "APEndOfLine", "APElement",
  "APElementWildCard", "APWildCard", "APOpenParenthesis",
  "APCloseParenthesis", "APOpenBracket", "APOpenAngleBracket",
  "APCloseBracket", "APCloseAngleBracket", "APOperatorAnd", "APOperatorOr",
  "APColon", "APAntechamberARLevel", "APAntechamberRingMembership",
  "APAntechamberNoRingMembership", "APBond", "APNumber", "APTag",
  "APError", "APNotBondedToPrevious", "APBondedToPrevious",
  "APAtomTypeDefinition", "APOperatorNot", "$accept", "antechamberMatch",
  "antechamberLine", "typeName", "tagOrElement", "residueList",
  "focusAtomicNumber", "focusNumberOfAttachedAtoms",
  "focusNumberOfAttachedHydrogenAtoms",
  "focusNumberOfElectronWithdrawingAtoms", "bracketedAtomicProperty",
  "chemicalEnvironment", "bondDefinitions", "logOp", "atomicTest",
  "branches", "branchesList", "chain", "chemicalEnvironmentAtom",
  "element", "bracketedAtomPropertyOrNull", "tagNameOrNull", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281
};
# endif

#define YYPACT_NINF -43

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-43)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
      -9,    10,    28,   -43,   -43,   -43,    50,    30,   -43,   -43,
     -43,    -1,     0,   -43,   -43,   -43,     3,    10,   -43,   -43,
     -43,   -43,     4,   -43,   -43,   -43,   -43,     5,   -43,   -43,
     -43,    24,   -43,   -43,    22,    40,   -43,    16,   -43,    11,
     -16,    37,   -43,   -43,   -43,    53,    -3,   -43,   -43,   -43,
     -43,    25,   -43,   -43,    22,    22,   -43,   -43,     9,   -43,
      42,     6,   -43,    44,    57,   -43,    47,   -43,   -43,    53,
     -43,    22,    54,    52,    43,   -43,   -43,    41,    52,    45,
     -43,    55,   -43,   -43,    56,    46,   -43,   -43
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,    17,    16,     2,     0,    12,     1,    11,
      18,     0,     0,    10,    19,    20,     0,    13,    15,     9,
      21,    22,     0,    14,     8,    23,    24,     0,     7,    25,
      26,     0,     6,    27,     0,     0,    39,    40,    43,    35,
       0,     0,    32,     5,    29,     0,     0,    30,    41,    37,
      38,     0,    36,    28,     0,     0,    51,    52,     0,    45,
      47,    53,     4,     0,     0,    42,    33,    34,    44,     0,
      48,     0,    53,    55,     0,     3,    46,     0,    55,     0,
      50,     0,    54,    49,     0,     0,    56,    31
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -43,   -43,   -43,   -43,    20,   -43,   -43,   -43,   -43,   -43,
     -43,   -43,   -43,   -42,   -43,     7,   -43,     2,   -43,   -43,
       1,    -6
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     2,     5,     6,     7,    11,    16,    22,    27,    31,
      35,    46,    64,    41,    42,    47,    58,    59,    60,    61,
      73,    80
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      62,    51,    13,    52,     3,    14,    19,    24,    28,    20,
      25,    29,    66,    67,     3,    71,     1,    68,    63,    15,
      17,     4,    69,    21,    26,    30,    72,    32,     8,    77,
      33,     4,    18,    34,    49,    50,    48,    23,    36,    37,
      38,    39,    40,    43,    12,    65,    44,    45,    53,    45,
      54,    55,    82,     9,    54,    55,    10,    56,    57,    74,
      75,    55,    79,    71,    81,    87,    84,    70,    86,     0,
      85,    76,    83,    78
};

static const yytype_int8 yycheck[] =
{
       3,    17,     3,    19,     4,     6,     3,     3,     3,     6,
       6,     6,    54,    55,     4,     9,    25,     8,    21,    20,
      20,    21,    13,    20,    20,    20,    20,     3,     0,    71,
       6,    21,    12,     9,    23,    24,    20,    17,    16,    17,
      18,    19,    20,     3,    14,    20,     6,     7,    11,     7,
      13,    14,    11,     3,    13,    14,     6,     4,     5,    15,
       3,    14,    10,     9,    21,    19,    21,    60,    12,    -1,
      15,    69,    78,    72
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    25,    28,     4,    21,    29,    30,    31,     0,     3,
       6,    32,    14,     3,     6,    20,    33,    20,    31,     3,
       6,    20,    34,    31,     3,     6,    20,    35,     3,     6,
      20,    36,     3,     6,     9,    37,    16,    17,    18,    19,
      20,    40,    41,     3,     6,     7,    38,    42,    20,    23,
      24,    17,    19,    11,    13,    14,     4,     5,    43,    44,
      45,    46,     3,    21,    39,    20,    40,    40,     8,    13,
      42,     9,    20,    47,    15,     3,    44,    40,    47,    10,
      48,    21,    11,    48,    21,    15,    12,    19
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    27,    28,    29,    29,    29,    29,    29,    29,    29,
      29,    29,    30,    30,    30,    30,    31,    31,    32,    33,
      33,    34,    34,    35,    35,    36,    36,    37,    37,    38,
      38,    39,    40,    40,    40,    41,    41,    41,    41,    41,
      41,    41,    41,    41,    42,    43,    43,    44,    44,    45,
      45,    46,    46,    47,    47,    48,    48
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,    10,     9,     8,     7,     6,     5,     4,
       3,     2,     1,     3,     4,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     1,
       1,     5,     1,     3,     3,     1,     2,     2,     2,     1,
       1,     2,     3,     1,     3,     1,     3,     1,     2,     4,
       3,     1,     1,     0,     3,     0,     3
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
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, gaff_SParserParam* data)
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, gaff_SParserParam* data)
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
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, gaff_SParserParam* data)
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, gaff_SParserParam* data)
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
          case 29: /* antechamberLine  */
#line 122 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).antechamberRoot)!=NULL) delete (((*yyvaluep).antechamberRoot)); }
#line 1166 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 32: /* residueList  */
#line 125 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).residueList)!=NULL) delete (((*yyvaluep).residueList)); }
#line 1172 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 37: /* bracketedAtomicProperty  */
#line 126 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).logical)!=NULL) delete (((*yyvaluep).logical)); }
#line 1178 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 38: /* chemicalEnvironment  */
#line 123 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).bondListMatchNode)!=NULL) delete (((*yyvaluep).bondListMatchNode)); }
#line 1184 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 39: /* bondDefinitions  */
#line 124 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).afterMatchBondTest)!=NULL) delete (((*yyvaluep).afterMatchBondTest)); }
#line 1190 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 40: /* logOp  */
#line 127 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).logical)!=NULL) delete (((*yyvaluep).logical)); }
#line 1196 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 41: /* atomicTest  */
#line 128 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).atomTest)!=NULL) delete (((*yyvaluep).atomTest)); }
#line 1202 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 42: /* branches  */
#line 129 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).bondListMatchNode)!=NULL) delete (((*yyvaluep).bondListMatchNode)); }
#line 1208 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 43: /* branchesList  */
#line 130 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).bondListMatchNode)!=NULL) delete (((*yyvaluep).bondListMatchNode)); }
#line 1214 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 44: /* chain  */
#line 131 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).bondListMatchNode)!=NULL) delete (((*yyvaluep).bondListMatchNode)); }
#line 1220 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 45: /* chemicalEnvironmentAtom  */
#line 132 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).antechamberBondTest)!=NULL) delete (((*yyvaluep).antechamberBondTest)); }
#line 1226 "_gaff_Parser.cpp" /* yacc.c:1257  */
        break;

    case 47: /* bracketedAtomPropertyOrNull  */
#line 133 "gaff_Parser.yy" /* yacc.c:1257  */
      { if (((*yyvaluep).logical)!=NULL) delete (((*yyvaluep).logical)); }
#line 1232 "_gaff_Parser.cpp" /* yacc.c:1257  */
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
yyparse (gaff_SParserParam* data)
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
#line 142 "gaff_Parser.yy" /* yacc.c:1661  */
    { data->expression = (yyvsp[0].antechamberRoot)->_obj; LOG(BF("Got AntechamberMatch tree=%0X") % (yyvsp[0].antechamberRoot) ); }
#line 1526 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 3:
#line 156 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( (yyvsp[-8].residueList)->_obj, (yyvsp[-7].ival), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].logical)->_obj );
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-9].carr),focus,(yyvsp[-2].bondListMatchNode)->_obj,(yyvsp[-1].afterMatchBondTest)->_obj)); 
                }
#line 1537 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 4:
#line 171 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( (yyvsp[-7].residueList)->_obj, (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].logical)->_obj );
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-8].carr),focus,(yyvsp[-1].bondListMatchNode)->_obj)); 
                }
#line 1548 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 5:
#line 185 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( (yyvsp[-6].residueList)->_obj, (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].logical)->_obj );
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-7].carr),focus)); 
                }
#line 1559 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 6:
#line 198 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
		    chem::AtomOrBondMatchNode_sp temp;
		    temp = chem::AtomOrBondMatchNode_O::_nil;
                    focus = chem::AntechamberFocusAtomMatch_O::create( (yyvsp[-5].residueList)->_obj, (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].ival), temp );
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-6].carr),focus)); 
                }
#line 1572 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 7:
#line 212 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( (yyvsp[-4].residueList)->_obj, (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].ival), -1);
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-5].carr),focus)); 
                }
#line 1583 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 8:
#line 223 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
		    focus = chem::AntechamberFocusAtomMatch_O::create( (yyvsp[-3].residueList)->_obj, (yyvsp[-2].ival), (yyvsp[-1].ival), -1, -1);
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-4].carr),focus));
                }
#line 1594 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 9:
#line 233 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( (yyvsp[-2].residueList)->_obj, (yyvsp[-1].ival), -1, -1, -1 );
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-3].carr),focus)); 
                }
#line 1605 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 10:
#line 242 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
                    focus = chem::AntechamberFocusAtomMatch_O::create( (yyvsp[-1].residueList)->_obj, -1, -1, -1, -1);
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-2].carr),focus)); 
                }
#line 1616 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 11:
#line 250 "gaff_Parser.yy" /* yacc.c:1661  */
    {
		    _BLOCK_TRACE("Got antechamber line");
                    chem::AntechamberFocusAtomMatch_sp focus;
		    focus = chem::AntechamberFocusAtomMatch_O::_nil;
                    (yyval.antechamberRoot) = new chem::Hold<chem::AntechamberRoot_O>(chem::AntechamberRoot_O::create((yyvsp[-1].carr),focus));
                }
#line 1627 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 12:
#line 259 "gaff_Parser.yy" /* yacc.c:1661  */
    { strcpy((yyval.carr),(yyvsp[0].carr)); }
#line 1633 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 13:
#line 260 "gaff_Parser.yy" /* yacc.c:1661  */
    { stringstream ss; ss << (yyvsp[-2].carr) << "." << (yyvsp[0].ival); strcpy((yyval.carr),ss.str().c_str()); }
#line 1639 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 14:
#line 261 "gaff_Parser.yy" /* yacc.c:1661  */
    { stringstream ss; ss << (yyvsp[-3].carr) << "." << (yyvsp[-1].ival) << (yyvsp[0].carr); strcpy((yyval.carr),ss.str().c_str()); }
#line 1645 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 15:
#line 262 "gaff_Parser.yy" /* yacc.c:1661  */
    { stringstream ss; ss << (yyvsp[-2].carr) << "." << (yyvsp[0].carr); strcpy((yyval.carr),ss.str().c_str()); }
#line 1651 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 16:
#line 266 "gaff_Parser.yy" /* yacc.c:1661  */
    { strcpy((yyval.carr),(yyvsp[0].carr)); }
#line 1657 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 17:
#line 267 "gaff_Parser.yy" /* yacc.c:1661  */
    { strcpy((yyval.carr),(yyvsp[0].carr)); }
#line 1663 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 18:
#line 274 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.residueList) = new chem::Hold<chem::ResidueList_O>(); }
#line 1669 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 19:
#line 279 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = -1; }
#line 1675 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 20:
#line 280 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 1681 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 21:
#line 284 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = -1; }
#line 1687 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 22:
#line 285 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 1693 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 23:
#line 289 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = -1; }
#line 1699 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 24:
#line 290 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 1705 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 25:
#line 294 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = -1; }
#line 1711 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 26:
#line 295 "gaff_Parser.yy" /* yacc.c:1661  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 1717 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 27:
#line 299 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.logical) = new chem::Hold<chem::Logical_O>(); }
#line 1723 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 28:
#line 300 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); 
		LOG(BF("Got bracketedAtomicProperty") ); // vp0(( "Got bracketedAtomicProperty" ));
		(yyval.logical) = new chem::Hold<chem::Logical_O>((yyvsp[-1].logical)->_obj); }
#line 1731 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 29:
#line 307 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.bondListMatchNode) = new chem::Hold<chem::BondListMatchNode_O>(); }
#line 1737 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 30:
#line 308 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.bondListMatchNode) = new chem::Hold<chem::BondListMatchNode_O>((yyvsp[0].bondListMatchNode)->_obj); }
#line 1743 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 31:
#line 314 "gaff_Parser.yy" /* yacc.c:1661  */
    {_BT(); (yyval.afterMatchBondTest) = new chem::Hold<chem::AfterMatchBondTest_O>(chem::AfterMatchBondTest_O::create( (yyvsp[-4].carr), (yyvsp[-2].carr), (yyvsp[0].benum) )); }
#line 1749 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 32:
#line 318 "gaff_Parser.yy" /* yacc.c:1661  */
    {_BT(); (yyval.logical) = new chem::Hold<chem::Logical_O>(chem::Logical_O::create( chem::logIdentity, (yyvsp[0].atomTest)->_obj )); }
#line 1755 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 33:
#line 319 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT();(yyval.logical) = new chem::Hold<chem::Logical_O>(chem::Logical_O::create(chem::logHighPrecedenceAnd, (yyvsp[-2].logical)->_obj, (yyvsp[0].logical)->_obj )); }
#line 1761 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 34:
#line 320 "gaff_Parser.yy" /* yacc.c:1661  */
    {_BT(); (yyval.logical) = new chem::Hold<chem::Logical_O>(chem::Logical_O::create(chem::logOr, (yyvsp[-2].logical)->_obj, (yyvsp[0].logical)->_obj )); }
#line 1767 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 35:
#line 323 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPInBond, (yyvsp[0].benum), 1)); }
#line 1773 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 36:
#line 324 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPInBond, (yyvsp[0].benum), (yyvsp[-1].ival))); }
#line 1779 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 37:
#line 325 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPNotBondedToPrevious,(yyvsp[-1].benum))); }
#line 1785 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 38:
#line 326 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPBondedToPrevious,(yyvsp[-1].benum))); }
#line 1791 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 39:
#line 327 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPArLevel, (yyvsp[0].ival) )); }
#line 1797 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 40:
#line 328 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPRingSize, 1 )); }
#line 1803 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 41:
#line 329 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPRingSize, (yyvsp[0].ival) )); }
#line 1809 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 42:
#line 330 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPRingSize, (yyvsp[0].ival), (yyvsp[-2].ival) )); }
#line 1815 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 43:
#line 331 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.atomTest) = new chem::Hold<chem::AtomTest_O>(chem::AtomTest_O::create(chem::SAPNoRing)); }
#line 1821 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 44:
#line 335 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.bondListMatchNode)=(yyvsp[-1].bondListMatchNode); }
#line 1827 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 45:
#line 338 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.bondListMatchNode) = new chem::Hold<chem::BondListMatchNode_O>((yyvsp[0].bondListMatchNode)->_obj); }
#line 1833 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 46:
#line 339 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.bondListMatchNode) = new chem::Hold<chem::BondListMatchNode_O>(chem::Branch_O::create((yyvsp[-2].bondListMatchNode)->_obj,(yyvsp[0].bondListMatchNode)->_obj)); }
#line 1839 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 47:
#line 343 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT();(yyval.bondListMatchNode) = new chem::Hold<chem::BondListMatchNode_O>(chem::Chain_O::create((yyvsp[0].antechamberBondTest)->_obj));}
#line 1845 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 48:
#line 344 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.bondListMatchNode) = new chem::Hold<chem::BondListMatchNode_O>(chem::Chain_O::create((yyvsp[-1].antechamberBondTest)->_obj,(yyvsp[0].bondListMatchNode)->_obj)); }
#line 1851 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 49:
#line 352 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); 
                    (yyval.antechamberBondTest) = new chem::Hold<chem::AntechamberBondTest_O>(chem::AntechamberBondTest_O::create((yyvsp[-3].carr),(yyvsp[-2].ival),(yyvsp[-1].logical)->_obj,(yyvsp[0].carr))); 
                }
#line 1859 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 50:
#line 357 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); 
                    (yyval.antechamberBondTest) = new chem::Hold<chem::AntechamberBondTest_O>(chem::AntechamberBondTest_O::create((yyvsp[-2].carr), -1, (yyvsp[-1].logical)->_obj, (yyvsp[0].carr)));
                }
#line 1867 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 51:
#line 363 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); strcpy((yyval.carr),(yyvsp[0].carr)); }
#line 1873 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 52:
#line 364 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); strcpy((yyval.carr),(yyvsp[0].carr)); }
#line 1879 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 53:
#line 367 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.logical) = new chem::Hold<chem::Logical_O>(); }
#line 1885 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 54:
#line 368 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); (yyval.logical) = new chem::Hold<chem::Logical_O>((yyvsp[-1].logical)->_obj); }
#line 1891 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 55:
#line 372 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); strcpy((yyval.carr),""); }
#line 1897 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;

  case 56:
#line 373 "gaff_Parser.yy" /* yacc.c:1661  */
    { _BT(); strcpy((yyval.carr),(yyvsp[-1].carr)); }
#line 1903 "_gaff_Parser.cpp" /* yacc.c:1661  */
    break;


#line 1907 "_gaff_Parser.cpp" /* yacc.c:1661  */
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
#line 380 "gaff_Parser.yy" /* yacc.c:1906  */


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
		    LEXPRINT("APOperatorAnd character");
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
{
    gaff_SParserParam p;
    stringstream sin(input);
    p.expression = _Nil<_O>();
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



    

