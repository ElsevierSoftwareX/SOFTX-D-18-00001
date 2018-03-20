%{

  /** @file ergo_input_processor.c Parses the input.
      Uses bison code generator to generate the parses.
  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ergo_scripted.h"
#define YYERROR_VERBOSE
int yylex(void);
int yyerror(const char *s);
static const char *last_token = NULL;

%}

%union {
  double num;     /* for returning numbers */
  char str[256];  /* for returning strings */
  struct variable *var; /* for returning lvalues */
}
  
%token <num> NUMBER
%token DOT
%token <str> SYMBOL EQUAL STRING EOFTAG GETEXC GETPOL K_ALL HELP MOLTAG GHOSTTAG MOLDAL QUIT RUNTAG SYSTEM GHOST ANGSTROM PRECISION RANGE WARRANTY
%token LIST_DFT_FUNCS
%token IS_CHT_USED
%token SET_NTHREADS
%token	PLUS	MINUS	TIMES	DIVIDE	POWER
%token	LEFT_PARENTHESIS	RIGHT_PARENTHESIS
%token	EOL

%type <num> Expression
%type <var> Lvalue

%left	PLUS	MINUS
%left	TIMES	DIVIDE
%left	NEG
%right	POWER

%start Input
%%

Input:
        |  Line
	|  Line EOL Input
	;

Line:
         /* Empty */
	| Assignment 
        | Command
        | error { if(!ergo_scanner_reading_stdin) { yyerror("Aborted."); YYABORT; } }
	;

Assignment:
          Lvalue EQUAL Expression { es_assign_num($1, $3);}
        | Lvalue EQUAL STRING     { es_assign_str($1, $3);}
        ;

Lvalue:
          SYMBOL            { $$=es_find_var(NULL, $1); 
                              if(!$$) { last_token = $1;
                                        yyerror("Unknown variable");YYERROR; }}
        | Lvalue DOT SYMBOL { $$=es_find_var($1,   $3);
                              if(!$$) { last_token = $3;
                                        yyerror("Unknown variable");YYERROR;}}
        ;

Command: 
          MOLTAG EOL Molinput { es_mol_commit(); }
        | MOLTAG ANGSTROM EOL Molinput { es_mol_commit(); }
        | GHOSTTAG EOL Molinput { es_mol_commit(); }
        | GHOSTTAG ANGSTROM EOL Molinput { es_mol_commit(); }
        | HELP                { es_print_help(); }
        | HELP Lvalue         { es_print_help_var($2); }
        | LIST_DFT_FUNCS      { es_print_list_dft_funcs(); }
        | IS_CHT_USED         { es_print_is_cht_used(); }
        | PRECISION           { es_print_precision(); }
        | MOLDAL STRING       { if(es_mol_read_molecule($2,MOL_MAIN)) { 
                              yyerror("Reading MOLECULE failed"); YYERROR;} }
        | MOLDAL GHOST STRING { if(es_mol_read_molecule($3, MOL_GHOST)) { 
                              yyerror("Reading GHOST MOLECULE failed"); YYERROR;} }
        | RUNTAG STRING     { if(es_run($2, 0)) {
                              yyerror("RUN failed"); YYERROR;} }
        | SET_NTHREADS LEFT_PARENTHESIS Expression RIGHT_PARENTHESIS {
          if(es_set_nthreads($3)) { yyerror("setNThreads failed"); YYERROR;} }
        | SET_NTHREADS LEFT_PARENTHESIS STRING RIGHT_PARENTHESIS {
          if(es_set_nthreads_string($3)) { yyerror("setNThreads failed"); YYERROR;} }
        | SYSTEM STRING { puts($2); 
                   if(system($2) != 0) {yyerror("system() failed"); YYERROR;} }
        | WARRANTY { es_warranty(); }
        | GETEXC STRING NUMBER { if(es_getexc($2, $3)) {
                              yyerror("get_excited_state failed"); YYERROR;} }
        | GETPOL STRING STRING NUMBER { if(es_get_polarisability($2, $3, $4)) {
                              yyerror("get_polarisability failed"); YYERROR;} }
        | GETPOL STRING K_ALL NUMBER { if(es_get_polarisability($2, NULL,$4)) {
                              yyerror("get_polarisability failed"); YYERROR;} }
        | QUIT { YYACCEPT; }
        | RANGE NUMBER EQUAL NUMBER NUMBER STRING {
                         if(!es_assign_range(MOL_MAIN,$2, $4, $5, $6)) {
                            yyerror("Invalid main basis set range");YYERROR;}
                         }
        | RANGE GHOST NUMBER EQUAL NUMBER NUMBER STRING {
                         if(!es_assign_range(MOL_GHOST,$3, $5, $6, $7)) {
                            yyerror("Invalid ghost basis set range");YYERROR;}
                         }
        ;

Molinput:
        EOFTAG
        | Molline Molinput
        ;

Molline:
        SYMBOL NUMBER NUMBER NUMBER EOL { es_add_atom($1, $2, $3, $4); }
        ;

Expression:
          NUMBER			{ $$=$1;  }
        | Expression PLUS Expression	{ $$=$1+$3; }
	| Expression MINUS Expression	{ $$=$1-$3; }
	| Expression TIMES Expression	{ $$=$1*$3; }
	| Expression DIVIDE Expression	{ $$=$1/$3; }
	| MINUS Expression %prec NEG	{ $$=-$2; }
	| Expression POWER Expression	{ $$=pow($1,$3); }
	| LEFT_PARENTHESIS Expression RIGHT_PARENTHESIS	{ $$=$2; }
	;

%%

YYSTYPE yylval;
int ergo_scanner_lineno = 1;
int ergo_scanner_reading_stdin = 0;

int yyerror(const char *s) {
  if (last_token) {
    printf("line %d: %s at '%s'\n",ergo_scanner_lineno, s, last_token);
    last_token = NULL;
  } else {
    printf("line %d: %s\n",ergo_scanner_lineno, s);
  }
  return !ergo_scanner_reading_stdin;
}

#ifdef SCANNER_TEST
int main(void) {
  yyparse();
}
#endif
