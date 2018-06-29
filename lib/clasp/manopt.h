#ifndef MANOPT_H
#define MANOPT_H
/*
 *
 *	manopt.h
 *  declarartions for the option manager    
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 09/01/2008 11:13:03 AM CEST  
 *
 *  Revision of last commit: 
 *  $Rev: 113 $
 *  $Author: steve $
 *  $Date: 2009-07-01 09:45:02 +0200 (Wed, 01 Jul 2009) $
 *
 *
 *  $Id: manopt.h 113 2009-07-01 07:45:02Z steve $ 
 *  $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/manopt.h $
 *
 */

#define MANOPT_MAXSYNOPSIS 10000

typedef enum{
  FLAG,
  REQSTRINGOPT,
  REQCHAROPT,
  REQINTOPT,
  REQUINTOPT,
  REQDBLOPT,
  MANOPT_ENUMREQUIRED,    /*marker to distinguish required and optional*/
  FILEOPT,
  STRINGOPT,
  CHAROPT,
  INTOPT,
  UINTOPT,
  DBLOPT,
  INTRANGEOPT,
  UINTRANGEOPT,
  DBLRANGEOPT,
  LISTOPT,
  SELECTOPT,
  MANOPT_ENUMSIZE,     /*end of enumeration*/
  MANOPT_BLOCKSEPARATOR
} manopt_type;


typedef struct { 
  char *flagname;
  int noofvalues;
  char **values;
} manopt_arg;

typedef struct {
  int noofargs;
  manopt_arg* args;
} manopt_argset;

typedef struct {
  char shortopt;
  char *longopt;
  char *argdesc;
  char *helpmsg;
  char *defaultval;
  unsigned char set;
  unsigned char required;
  manopt_type type;
  void *constraint; 
  manopt_arg arg;
  void *reg_var;
} manopt_option;

typedef struct {
  char *call;
  char *unflagged;
  char *references;
  char *bugs;
  char *version;
  char *description;
  int noofopts;
  manopt_option *opts;
  void *mutually_exclusive_opts;
} manopt_optionset; 

typedef struct {
  int maxlength;
  int minlength;
  int noofitems;
  char** items;
} manopt_listconstraint;

typedef struct {
  int max;
  int min;
  int diff;
} manopt_intconstraint;

typedef struct {
  unsigned int max;
  unsigned int min;
  unsigned int diff;
} manopt_uintcontstraint;

typedef struct {
  double max;
  double min;
  double diff;
} manopt_dblconstraint;

int
manopt_parse_commandline(manopt_argset* argset, 
    int argc, 
    char **argv);

void 
manopt(manopt_optionset* set, 
    manopt_type type, 
    unsigned char required,
    char shortopt, 
    char *longopt, 
    char *helpmsg, 
    char *argdesc,
    void *constraints,
    void *reg_var);

manopt_arg*
manopt_getopts(manopt_optionset* set, 
    int argc, 
    char **argv);

void
manopt_destructarg(manopt_arg *arg);

void
manopt_destructoptionset(manopt_optionset *set);
 
void
manopt_dumpoptionset(manopt_optionset *set);

void
manopt_helpmsg(manopt_optionset *set);
 
void
manopt_help(manopt_optionset *set, const char *fmt, ...);

void
manopt_initoptionset(manopt_optionset *set, 
    char *call, char *unflagged, char *description, char *references, char *version, char *bugs);

manopt_option*
manopt_longopt(manopt_optionset *set, char *longopt);

manopt_option*
manopt_shortopt(manopt_optionset *set, char shortopt);

unsigned char
manopt_isset(manopt_optionset *set, char shortopt, char *longopt);
 
void
manopt_blockseparator(manopt_optionset *set, char *blockname);

manopt_arg*
manopt_getarg(manopt_optionset *set, char shortopt, char *longopt);

unsigned char isfloat(char *s);
unsigned char isint(char *s);
#endif

