%option noyywrap
%option nounput
%option noinput
%{

/** @file ergo_input_parser.c The lex input parser. Defines the grammar */
/* for ergo input files. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "ergo_input_processor.h"
#include "ergo_scripted.h"
#define YY_DECL int yylex(void)
YY_DECL;

/* Some flex program generate external functions that are not used by us.
   We declare apriopriate prototypes to hush the warnings...
*/
/* ELIAS NOTE 2013-06-12: In Fedora 18 there were problems with compiler errors like "conflicting types for ‘yyget_leng’".
   To fix that, the declarations in the nine following lines with declarations of yyget_leng, yyget_in etc were commented out.
   Probably those lines should not have been needed anyway, judging by the above old comment saying that the lines were added "to hush the warnings...".
*/
/*
static int yyget_leng(void);
static FILE *yyget_in  (void);
static FILE *yyget_out (void);
static char *yyget_text(void);
static void yyset_in (FILE *  in_str);
static void yyset_out(FILE *  out_str);
static int yyget_debug(void);
static void yyset_debug(int  bdebug);
static int yylex_destroy(void);
*/
extern int fileno(FILE *stream);
/* Following needed for flex-2.5.33-5 as in F7. */
void yyset_lineno (int  line_number );
int yyget_lineno  (void);
%}
/* parse the input and divide it into basic building blocks: tokens */

white		[ \t]+
comment1         ^#.*
comment2         ^;.*

digit		[0-9]
integer 	{digit}+
exponent	[eE][+-]?{integer}

real		-?{integer}("."{integer})?{exponent}?
runtag          "run"
moltag          "molecule_inline"
ghosttag        "ghost_inline"
moldal          "molecule"
getexc          "get_excited_state"
getpol          "get_polarisability"
system          "system"
warranty        "warranty"
symbol          [A-Za-z_][0-9A-Za-z_]*
string          \"[^"]*\"

%%
{white} 	{ /* We ignore white characters. */ }
{comment1} 	{ /* We ignore comments, too. */  }
{comment2} 	{ /* We ignore comments, too. */  }

{real}		{ yylval.num =atof(yytext); return(NUMBER); }
"."             return(DOT);
"+"		return(PLUS);
"-"		return(MINUS);
^EOF$           return EOFTAG;

"*"		return(TIMES);
"/"		return(DIVIDE);

"^"		return(POWER);

"("		return(LEFT_PARENTHESIS);
")"		return(RIGHT_PARENTHESIS);
"="             return(EQUAL);
"\n"     	{ if(ergo_scanner_reading_stdin) printf("> ");
                  ergo_scanner_lineno++;return(EOL); }
";"             return EOL;
"ghost"         return GHOST;
"help"          return HELP;
"list_dft_funcs" return LIST_DFT_FUNCS;
"is_cht_used"    return IS_CHT_USED;
"precision"     return PRECISION;
"range"         return RANGE;
"Angstrom"      { es_mol_unit_angstrom(); return ANGSTROM; }
"set_nthreads"  return SET_NTHREADS;
{getexc}          return GETEXC;
{getpol}          return GETPOL;
"all"             return K_ALL;
{moltag}        { es_mol_begin(MOL_MAIN);  return MOLTAG; }
{ghosttag}      { es_mol_begin(MOL_GHOST); return MOLTAG; }
{moldal}          return MOLDAL;
"quit"            return QUIT;
{runtag}          return RUNTAG;
{system}          return SYSTEM;
{warranty}        return WARRANTY;


{string}        { unsigned l = strlen(yytext)-2; /* skip quotes */
                  if(l>=sizeof(yylval.str))l=sizeof(yylval.str)-1;
                  strncpy(yylval.str, yytext+1, l); yylval.str[l] = '\0';
                  return(STRING); }
{symbol} { strncpy(yylval.str, yytext, sizeof(yylval.str)-1);
 yylval.str[sizeof(yylval.str)-1] = '\0';
 return (SYMBOL); }
