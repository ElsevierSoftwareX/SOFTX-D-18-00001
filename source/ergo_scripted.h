/* Ergo, version 3.6, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2017 Elias Rudberg, Emanuel H. Rubensson, Pawel Salek,
 * and Anastasia Kruchinina.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

#ifndef ERGO_SCRIPTED_HEADER
#define ERGO_SCRIPTED_HEADER
/** @file ergo_scripted.h contains ergo input processor definitions.
 * Currently it requires access to all the modules. It is trivial to
 * make all the modules to register modifable variables etc
 * settings. It is more complex to expose access to routines because
 * some stub routines placing arguments on stack may be needed. Unless
 * we skip arguments and always require that modules must be set up
 * and then a "run" method is called.
 * 
 */
#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

extern int ergo_scanner_lineno;
extern int ergo_scanner_reading_stdin;

/** VarType defines recognized variable types. */
enum VarType { VAR_STRING, VAR_FLOAT, VAR_INT, VAR_LIST };

/** describes a variable recognized by the scripting system. */
struct variable {
  const char *name;
  const char *description;
  union {
    char   *str;
    double    num;
    int       vint;  
    struct variable *list;
  } v;
  enum VarType type;
  struct variable *next;
};

/** MolType decides whether molecule data access routines modify the
    main or the ghost molecule. */
enum MolType { MOL_MAIN, MOL_GHOST };

EXTERN_C void es_assign_num(struct variable *var, double val);
EXTERN_C void es_assign_str(struct variable *var, const char*str);
EXTERN_C int es_assign_range(enum MolType mt, int rangeNo,
                             int start, int cnt, const char *name);
EXTERN_C struct variable *es_find_var(struct variable *root, const char *name);

EXTERN_C void es_mol_begin(enum MolType moleculeClass);
EXTERN_C void es_add_atom(const char *name, double x, double y, double z);
EXTERN_C void es_mol_commit(void);
EXTERN_C void es_mol_unit_angstrom(void);

EXTERN_C int es_mol_read_molecule(const char *fname,
                                  enum MolType moleculeClass);

EXTERN_C void es_print_help();
EXTERN_C void es_print_help_var(const struct variable *root);
EXTERN_C void es_print_list_dft_funcs();
EXTERN_C void es_print_is_cht_used();
EXTERN_C void es_print_precision();
EXTERN_C int es_run(const char *mode, int save_pot);
EXTERN_C void es_warranty(void);
EXTERN_C int es_getexc(const char *mode, int modes);
EXTERN_C int es_get_polarisability(const char *mode, const char *opname,
				   double freq);

EXTERN_C int es_set_nthreads(int nThreads);
EXTERN_C int es_set_nthreads_string(const char *str);
extern int es_quit;

#endif /* ERGO_SCRIPTED_HEADER */
