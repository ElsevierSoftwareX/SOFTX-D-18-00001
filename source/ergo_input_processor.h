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

#ifndef YY_YY_ERGO_INPUT_PROCESSOR_H_INCLUDED
# define YY_YY_ERGO_INPUT_PROCESSOR_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NUMBER = 258,
    DOT = 259,
    SYMBOL = 260,
    EQUAL = 261,
    STRING = 262,
    EOFTAG = 263,
    GETEXC = 264,
    GETPOL = 265,
    K_ALL = 266,
    HELP = 267,
    MOLTAG = 268,
    GHOSTTAG = 269,
    MOLDAL = 270,
    QUIT = 271,
    RUNTAG = 272,
    SYSTEM = 273,
    GHOST = 274,
    ANGSTROM = 275,
    PRECISION = 276,
    RANGE = 277,
    WARRANTY = 278,
    LIST_DFT_FUNCS = 279,
    IS_CHT_USED = 280,
    SET_NTHREADS = 281,
    PLUS = 282,
    MINUS = 283,
    TIMES = 284,
    DIVIDE = 285,
    POWER = 286,
    LEFT_PARENTHESIS = 287,
    RIGHT_PARENTHESIS = 288,
    EOL = 289,
    NEG = 290
  };
#endif
/* Tokens.  */
#define NUMBER 258
#define DOT 259
#define SYMBOL 260
#define EQUAL 261
#define STRING 262
#define EOFTAG 263
#define GETEXC 264
#define GETPOL 265
#define K_ALL 266
#define HELP 267
#define MOLTAG 268
#define GHOSTTAG 269
#define MOLDAL 270
#define QUIT 271
#define RUNTAG 272
#define SYSTEM 273
#define GHOST 274
#define ANGSTROM 275
#define PRECISION 276
#define RANGE 277
#define WARRANTY 278
#define LIST_DFT_FUNCS 279
#define IS_CHT_USED 280
#define SET_NTHREADS 281
#define PLUS 282
#define MINUS 283
#define TIMES 284
#define DIVIDE 285
#define POWER 286
#define LEFT_PARENTHESIS 287
#define RIGHT_PARENTHESIS 288
#define EOL 289
#define NEG 290

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 18 "ergo_input_processor.y" /* yacc.c:1909  */

  double num;     /* for returning numbers */
  char str[256];  /* for returning strings */
  struct variable *var; /* for returning lvalues */

#line 130 "ergo_input_processor.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_ERGO_INPUT_PROCESSOR_H_INCLUDED  */
