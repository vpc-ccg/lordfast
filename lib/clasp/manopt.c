
/*
 *  manopt.c
 *  implementations for the option manager 
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 09/01/2008 11:12:56 AM CEST
 *
 *  Revision of last commit: 
 *  $Rev: 74 $
 *  $Author: steve $
 *  $Date: 2008-10-29 15:03:04 +0100 (Wed, 29 Oct 2008) $
 *
 *
 *  $Id: manopt.c 74 2008-10-29 14:03:04Z steve $ 
 *  $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/manopt.c $
 *
 */

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <ctype.h>
#include "manopt.h"
#include <sys/ioctl.h>


int
detectTerminalwidth(void) {
  struct winsize termwinsize;
  ioctl(0,TIOCGWINSZ, &termwinsize);
  return termwinsize.ws_col;
}

unsigned char
isfloat(char *s) {
  int i=0;
  int len=strlen(s);
  unsigned char dpt = 0;

  i += (s[i] == 43 || s[i] == 45) ? 1 : 0;
  while((s[i] >= 48 && s[i] <= 57) || (dpt=(!dpt && s[i]==46))) i++;

  return (len==i);
}

unsigned char
isint(char *s) {
  int i=0;
  int len=strlen(s);
  i += (s[i] == 43 || s[i] == 45) ? 1 : 0;
  while((s[i] >= 48 && s[i] <= 57)) i++;
  return (len==i);
}

void
manopt_usage(manopt_optionset *set) {
  unsigned int i=0,j=0,k,l,
      aptr = 0,
      msglen = 0,
      restlen = 0,
      lastspace = 0,
      offset =0,
      synopsislen=0,
      calllen=0;
  char shortopt = 0;
  char *longopt = NULL;
  unsigned char flags = 0;
  unsigned int maxarglen = 0;
  
  char string[2];
  char *synopsis;
  char **arg;
  char **msg;
  char *fill=NULL;
  char *call;
  int width = detectTerminalwidth();
  
  width = (width < 70) ? 70 : width;
  synopsis = malloc(sizeof(char)*MANOPT_MAXSYNOPSIS);
  call = malloc(sizeof(char)*MANOPT_MAXSYNOPSIS);
  call[0] = '\0',
  synopsis[0] = '\0' ;
  
  arg = malloc(sizeof(char*)*set->noofopts);
  msg = malloc(sizeof(char*)*set->noofopts);
  
  for(i=0;i<set->noofopts;i++) {
    arg[i] = malloc(sizeof(char)*MANOPT_MAXSYNOPSIS);
    arg[i][0] = 0;
    msg[i] = malloc(sizeof(char)*MANOPT_MAXSYNOPSIS);
    msg[i][0] = 0;
  }

  strcat(call, "usage: ");
  strcat(call, set->call);
  strcat(call, " ");
  calllen = strlen(call);
  if (calllen > 40) {
    strcat(synopsis, "\n");
    calllen = 20;
  }

  for(i=0; i < set->noofopts; i++) {
    shortopt = set->opts[i].shortopt;
    longopt = set->opts[i].longopt;
    if (set->opts[i].type == FLAG && shortopt) {
      if (!flags) {
        strcat(synopsis, "[-");
        flags = 1;
      }
      string[0] = set->opts[i].shortopt;
      string[1] = 0;
      strcat(synopsis, string);
   }
  }

  if(flags) {
    strcat(synopsis, "]\t");
  }

  flags = 0;
  for(i=0; i < set->noofopts; i++) {
    shortopt = set->opts[i].shortopt;
    longopt = set->opts[i].longopt;
    if (set->opts[i].type == MANOPT_BLOCKSEPARATOR) {
      strcat(arg[aptr], " [");
      strcat(arg[aptr], set->opts[i].longopt); 
      strcat(arg[aptr], "]");
      msg[aptr][0]=0;
      aptr++;
    }
    else if (set->opts[i].type != FLAG || 
        (set->opts[i].type == FLAG && !shortopt)) {
      if (!set->opts[i].required) {
        strcat(synopsis, "[");
        flags = 1;
      }
      if(shortopt) {
        strcat(arg[aptr], " ");
        strcat(synopsis, "-");
        strcat(arg[aptr], "-");
        string[0] = set->opts[i].shortopt;
        string[1] = 0;
        strcat(synopsis, string);
        strcat(arg[aptr],string);
        if(longopt) {
          strcat(arg[aptr], ",");
        }
      }
      if(longopt) {
        strcat(arg[aptr]," --");
        strcat(arg[aptr], longopt);
        if(!shortopt) {
          strcat(synopsis,"--");
          strcat(synopsis, longopt);
        }
      }

      if(set->opts[i].argdesc) {
        strcat(arg[aptr], " ");
        strcat(arg[aptr], set->opts[i].argdesc);
        strcat(synopsis, " ");
        strcat(synopsis, set->opts[i].argdesc); 
      }
      strcat(arg[aptr], " ");

      strcat(msg[aptr], set->opts[i].helpmsg);
      if (set->opts[i].defaultval) {
        strcat(msg[aptr], " (default:");
        strcat(msg[aptr], set->opts[i].defaultval); 
        strcat(msg[aptr], ")");
      }
      aptr++;

      if (!set->opts[i].required) {
        strcat(synopsis, "]\t");
        flags = 0;
      } else {
        strcat(synopsis, "\t");
      }
    } else { 
      string[0] = set->opts[i].shortopt;
      string[1] = 0;
      strcat(arg[aptr], " ");
      strcat(arg[aptr], "-");
      strcat(arg[aptr], string);
      if(longopt) {
        strcat(arg[aptr], ", ");
        strcat(arg[aptr], "--");
        strcat(arg[aptr], longopt);
      }
      strcat(arg[aptr], " ");
      strcat(msg[aptr], set->opts[i].helpmsg);
      aptr++;
    }
  }

  if(set->unflagged) {
    strcat(synopsis, set->unflagged);
    strcat(synopsis, "\t");
  }

  fill = realloc(fill, calllen*sizeof(char));
  for(i=0; i < calllen; i++) {
    fill[i]=' ';
  }

  synopsislen = strlen(synopsis);
  if(calllen+synopsislen > width) {
    l = (synopsislen)/(width-calllen)+1;
    offset =0;
    for(k=0; k < l; k++) {
      for(j=(k*(width-calllen)); 
          j < ((k+1)*(width-calllen))-1-offset && j < strlen(synopsis); j++) {
        if(synopsis[j]=='\t') {
          lastspace = j;
        }
      }  
      offset = (k+1)*(width-calllen)-lastspace;
      restlen = strlen(&synopsis[lastspace+1]);
      memmove(&synopsis[lastspace+2], &synopsis[lastspace+1],
          (restlen)*sizeof(char));
      synopsis[lastspace+2+restlen] = 0;
      synopsis[lastspace+1] = '\n';
    }
    for(j=0; j < strlen(synopsis);j++) {
      if(synopsis[j] == '\t') synopsis[j] = ' '; 
      if(synopsis[j] == '\n') {
        restlen = strlen(&synopsis[j+1]);
        memmove(&synopsis[j+1+calllen], &synopsis[j+1], (restlen)*sizeof(char));
        synopsis[j+1+calllen+restlen] = 0;
        memmove(&synopsis[j+1], fill, calllen*sizeof(char));
      }
    }
  } else {
  
    for(k=0; k < strlen(synopsis); k++) {
      if(synopsis[k] == '\t') synopsis[k] = ' '; 
    }
  }

  for(i=0; i < set->noofopts;  i++) { 
    maxarglen = strlen(arg[i]) > maxarglen ? strlen(arg[i]) : maxarglen;
  }
  maxarglen++;
  assert(maxarglen < 60);
  
  fill = realloc(fill, maxarglen*sizeof(char));
  for(i=0; i < maxarglen; i++) {
    fill[i]=' ';
  }

  for(i=0; i < set->noofopts; i++) {
    if((msglen=strlen(msg[i])) > width-maxarglen) {
      l = (msglen)/(width-maxarglen)+1;
      offset =0;
      for(k=0; k < l; k++) {
        for(j=(k*(width-maxarglen)); 
            j < ((k+1)*(width-maxarglen))-1-offset && j < strlen(msg[i])
            ; j++) {
          if(isspace(msg[i][j])) {
            lastspace = j; 
          }
        }
        if (j >= (k+1)*(width-maxarglen)-1-offset) {
        offset = (k+1)*(width-maxarglen)-lastspace;
        restlen = strlen(&msg[i][lastspace+1]);
        memmove(&msg[i][lastspace+2],&msg[i][lastspace+1],
            restlen*sizeof(char));
        msg[i][lastspace+2+restlen] = 0;
        msg[i][lastspace+1] = '\n';
        }
      }
      for(j=0; j < strlen(msg[i]); j++) {
        if(msg[i][j] == '\n') {
         restlen = strlen(&msg[i][j+1]);
         memmove(&msg[i][j+1+maxarglen],&msg[i][j+1], restlen*sizeof(char));
         msg[i][j+1+maxarglen+restlen] = 0;
         memmove(&msg[i][j+1], fill, maxarglen*sizeof(char));
        }
      }
    }
  }

  fprintf(stderr, "%s", call);
  fprintf(stderr, "%s\n" ,synopsis);
  fprintf(stderr, "%s\n", set->description);
  for(i=0; i < set->noofopts;  i++) { 
    fprintf(stderr, "%s", arg[i]);
    for(j=0; j < maxarglen-strlen(arg[i]); j++) {
      fprintf(stderr, " ");
    }
    fprintf(stderr, "%s\n", msg[i]);
  }
  fprintf(stderr, " [VERSION]\n%s\n", set->version);
  fprintf(stderr, " [BUGS]\n%s\n", set->bugs);
  fprintf(stderr, " [REFERENCES]\n%s\n", set->references);
  
  
  for(i=0;i<set->noofopts;i++) {
    free(arg[i]); 
    free(msg[i]); 
  }

  free(fill);
  free(call);
  free(synopsis);
  free(arg);
  free(msg);

}

void
manopt_help(manopt_optionset *set, const char *fmt, ...) {
  int ret;
  va_list ap;
  va_start(ap, fmt);

  fprintf(stderr, "%s: ", set->call);
  ret = vfprintf(stderr, fmt, ap);
  va_end(ap);
  manopt_usage(set);

  exit(-1);
  return;
}

void
manopt_initoptionset(manopt_optionset *set, 
    char *call, char *unflagged, 
    char *description,
    char *references, 
    char *version, 
    char *bugs) {
  
  set->call = call;
  set->description = description;
  set->unflagged = unflagged;
  set->references = references;
  set->bugs = bugs;
  set->version = version;
  set->opts = NULL;
  set->noofopts = 0;
  return;
}
void
manopt_initarg(manopt_arg* arg) {
  arg->flagname = NULL;
  arg->noofvalues = 0;
  arg->values = NULL;
  return;
}

void
manopt_initoption(manopt_option *opt) {
  opt->shortopt = 0;
  opt->longopt = NULL;
  opt->helpmsg = NULL;
  opt->type = 0;
  opt->required = 0;
  opt->constraint = NULL;
  opt->set = (unsigned char) 0;
  opt->defaultval = NULL;
  opt->reg_var = NULL;
  manopt_initarg(&(opt->arg));

  return;
}


void
manopt_destructarg(manopt_arg *arg) {
  free(arg->values);
  return;
}

void
manopt_destructoptionset(manopt_optionset *set) {
  int i;

  for(i=0; i < set->noofopts; i++) {
    if (set->opts[i].arg.noofvalues) {
      manopt_destructarg(&set->opts[i].arg);
    }
    if(set->opts[i].defaultval) {
        free(set->opts[i].defaultval);
        set->opts[i].defaultval = NULL;
    }
  }
  if (set->noofopts > 0) free(set->opts);
  return;
}

int
manopt_parse_commandline(manopt_argset* argset, int argc, char **argv) {
  int i, 
      cnt=0, 
      len=0, 
      offset=0;
  manopt_arg *arg=NULL;

  for (i=0; i < argc; i++) {
    /*if a number follows '-' expression is considered argument*/
    if(argv[i][0] == '-' && (argv[i][1] < 48 || argv[i][1] > 57)) {
      offset = (argv[i][1] == '-') ? 1 : 0;
      arg = realloc(arg, sizeof(manopt_arg)*(cnt+1));
      manopt_initarg(&arg[cnt]);
      len = strlen(&argv[i][offset+1])+1;
      if (len <= 0) {
        fprintf(stderr, "flaglen <= 0!");
        return 0;
      }
      arg[cnt].flagname=&argv[i][offset+1];
      cnt++;
    } else {
      if(cnt == 0) {
        arg = realloc(arg, sizeof(manopt_arg)*(cnt+1));
        manopt_initarg(&arg[cnt]); 
        cnt++;
      }       
      arg[cnt-1].values = realloc(arg[cnt-1].values, 
          sizeof(char*)*(arg[cnt-1].noofvalues+1)); 
      arg[cnt-1].values[arg[cnt-1].noofvalues] = (char*) argv[i];
      arg[cnt-1].noofvalues++;
    }
  }
  argset->noofargs = cnt;
  argset->args= arg;

  return 1;
}

void
manopt_blockseparator(manopt_optionset *set, char *blockname) {

  set->opts = realloc(set->opts, sizeof(manopt_option)*(set->noofopts+1));
  manopt_initoption(&set->opts[set->noofopts]);
  set->opts[set->noofopts].longopt = blockname;
  set->opts[set->noofopts].type = MANOPT_BLOCKSEPARATOR;
  set->noofopts++;
}

void
manopt(manopt_optionset *set,
    manopt_type type,
    unsigned char required,
    char shortopt, 
    char *longopt,
    char *helpmsg,
    char *argdesc,
    void *constraint,
    void *reg_var) {
  unsigned int *uintval, 
               *uintrangeval;
  char *charval, 
       **ptr;
  int i, 
      *intval, 
      *intrangeval;
  double *dblval, 
         *dblrangeval;

  for(i=0; i < set->noofopts; i++) {
    if (shortopt && shortopt == set->opts[i].shortopt) {
      fprintf(stderr, "shortopt %c already defined", shortopt);
      exit(-1);
    }
    if (longopt && !strcmp(longopt,set->opts[i].longopt)) {
      fprintf(stderr, "longopt %s already defined", longopt);
      exit(-1);
    }
  }

  set->opts = realloc(set->opts, sizeof(manopt_option)*(set->noofopts+1));
  manopt_initoption(&set->opts[set->noofopts]);
  set->opts[set->noofopts].argdesc = argdesc;
  set->opts[set->noofopts].shortopt = shortopt;
  set->opts[set->noofopts].longopt = longopt;
  set->opts[set->noofopts].helpmsg = helpmsg;
  set->opts[set->noofopts].type = type;
  set->opts[set->noofopts].required = required;
  set->opts[set->noofopts].constraint = constraint;
  set->opts[set->noofopts].reg_var = reg_var;
  set->noofopts++;

  if(reg_var) {
    set->opts[i].defaultval = malloc(sizeof(char)*MANOPT_MAXSYNOPSIS);
    set->opts[i].defaultval[0]=0;
    switch(type) {
      case CHAROPT:
        charval = (char*) set->opts[i].reg_var;
        sprintf(set->opts[i].defaultval, "%c", charval[0]); 
        break;
      case REQUINTOPT:
      case REQINTOPT:
        uintval = (unsigned int*) set->opts[i].reg_var;
        sprintf(set->opts[i].defaultval, "%d", uintval[0]);
        break;
      case UINTOPT:
      case INTOPT:
        intval = (int*) set->opts[i].reg_var;
        sprintf(set->opts[i].defaultval, "%d", intval[0]);
        break;
      case REQDBLOPT:
      case DBLOPT:
        dblval = (double*) set->opts[i].reg_var; 
        sprintf(set->opts[i].defaultval, "%f", dblval[0]);                 
        break;
      case REQSTRINGOPT:
      case STRINGOPT:
        ptr = (char**) set->opts[i].reg_var;
        if(ptr[0])  {
          sprintf(set->opts[i].defaultval, "\"%s\"", ptr[0]);                 
        } else {
          sprintf(set->opts[i].defaultval, "none");                 
        }
        break;
      case INTRANGEOPT:
        intrangeval = (int*) set->opts[i].reg_var;
        if (intrangeval) {
          sprintf(set->opts[i].defaultval, "[%d,%d]", 
              intrangeval[0], intrangeval[1]);
        }
        break;
      case UINTRANGEOPT:
        uintrangeval = (unsigned int*) set->opts[i].reg_var;
        if (uintrangeval) {
          sprintf(set->opts[i].defaultval, "[%d,%d]", 
              uintrangeval[0], uintrangeval[1]);
        }                
        break;
      case DBLRANGEOPT:
        dblrangeval = (double*) set->opts[i].reg_var;
        if (dblrangeval) {
          sprintf(set->opts[i].defaultval, "[%f,%f]", 
              dblrangeval[0], dblrangeval[1]);
        }
        break;
      default:
        free(set->opts[i].defaultval);
        set->opts[i].defaultval = NULL;
        break;
    }
  }

  return;
}

void
manopt_unflag(manopt_argset *argset, 
    int arg, 
    int offset) {
  int size = argset->args[arg].noofvalues-offset;

  argset->args[0].values = realloc(argset->args[0].values, 
      (argset->args[0].noofvalues+size)*sizeof(char*));
  memmove(&argset->args[0].values[argset->args[0].noofvalues],
      &argset->args[arg].values[offset], size*sizeof(char*));
  argset->args[0].noofvalues += size;
  argset->args[arg].values = realloc(argset->args[arg].values, 
      size*sizeof(char*));
  argset->args[arg].noofvalues = offset;

  return;
}


unsigned char
manopt_checkconstraint(manopt_optionset* optset, 
    int opt, 
    manopt_argset *argset, 
    int arg) {

  unsigned char lastarg = (unsigned char) (arg == argset->noofargs-1);   
  int noofvalues = argset->args[arg].noofvalues;
  int i, 
      j, 
      rintval = 0, 
      lintval = 0;
  void *constraint = optset->opts[opt].constraint;
  double ldblval = .0, 
         rdblval = .0;
  unsigned char valid_select = 0;

  manopt_dblconstraint *dblconstraint = NULL;
  manopt_intconstraint *intconstraint = NULL;
  manopt_listconstraint *listconstraint = NULL;

  switch(optset->opts[opt].type) {
    case FLAG:
      if(noofvalues > 0) {
        if(!lastarg) {
          manopt_help(optset, "flag %c (%s) with argument given\n", 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else {
          manopt_unflag(argset, arg, 0);
        }
      }
      break;
    case REQCHAROPT:
      if(noofvalues < 1) {
        manopt_help(optset, "option %c (%s) without required argument\n",
            optset->opts[opt].shortopt, optset->opts[opt].longopt);
      }
    case CHAROPT:
      if(noofvalues > 1) {
        if(!lastarg) {
          manopt_help(optset, "option with multiple arguments\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else {
          manopt_unflag(argset, arg, 1);
        }
      } else if(strlen(argset->args[arg].values[0]) > 1) {
        manopt_help(optset, "a char for option %c (%s) argument required\n",
            optset->opts[opt].shortopt, optset->opts[opt].longopt);
      }
      break;
    case REQSTRINGOPT:
      if (noofvalues < 1) {
        manopt_help(optset, "option %c (%s) without required argument\n",
            optset->opts[opt].shortopt, optset->opts[opt].longopt);
      }
    case FILEOPT:
    case STRINGOPT:
      if (noofvalues > 1) {
        if(!lastarg) {
          manopt_help(optset,"option %c (%s) with multiple arguments\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else {
          manopt_unflag(argset, arg, 1);
        }
      }
      break;
    case REQDBLOPT:
      if (noofvalues < 1) {
        manopt_help(optset, "option %c (%s) without required argument\n",
            optset->opts[opt].shortopt, optset->opts[opt].longopt);
      }
    case DBLOPT:
      if (noofvalues > 1) {
        if(!lastarg) {
          manopt_help(optset, "option %c (%s) with multiple arguments\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else {  
          manopt_unflag(argset, arg, 1);
        }
      } else if (noofvalues) {
        if (!isfloat(argset->args[arg].values[0]) || 
            (ldblval=atof(argset->args[arg].values[0])) == HUGE_VAL) {
          manopt_help(optset, "double '%s' argument for option %c (%s) out of range\n",
              argset->args[arg].values[0], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else {
          dblconstraint = (manopt_dblconstraint*) constraint;
          if (dblconstraint && 
              (ldblval > dblconstraint->max || ldblval < dblconstraint->min)) {
            manopt_help(optset, "double '%s' argument for option %c (%s) out of bounds\n",
              argset->args[arg].values[0], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
          }
        }
      }
      break;
    case REQINTOPT:
      if (noofvalues < 1) {
        manopt_help(optset, "option %c (%s) without required argument\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
      }
    case INTOPT:
      if (noofvalues > 1) {
        if(!lastarg) {
          manopt_help(optset, "option %c (%s) with multiple arguments\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else { 
          manopt_unflag(argset, arg, 1);
        }
      } else if (noofvalues) {
        if (!isint(argset->args[arg].values[0]) || 
            (lintval=atoi(argset->args[arg].values[0])) == INT_MIN || 
            lintval == INT_MAX) {
          manopt_help(optset, "int argument '%s' for option %c (%s) out of range\n",
              argset->args[arg].values[0],
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else { 
          intconstraint = (manopt_intconstraint*) constraint;
          if (intconstraint && 
              (lintval > intconstraint->max || lintval < intconstraint->min)) {
            manopt_help(optset, "int argument '%s' for option %c (%s) out of bounds\n",
              argset->args[arg].values[0],
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
          }           
        }       
      }
      break;
    case REQUINTOPT:
      if (noofvalues < 1) {
        manopt_help(optset, "option %c (%s) without required argument\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
      }
    case UINTOPT:
      if (noofvalues > 1) {
        if(!lastarg) {
          manopt_help(optset, "option %c (%s) with multiple arguments\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else { 
          manopt_unflag(argset, arg, 1);
        }
      } else if (noofvalues) {
        if (!isint(argset->args[arg].values[0]) || 
            (lintval=atoi(argset->args[arg].values[0])) < 0 || 
            lintval == INT_MAX) {
          manopt_help(optset, "unsigned int argument '%s' for option %c (%s) out of range\n",
              argset->args[arg].values[0], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);

        } else { 
          intconstraint = (manopt_intconstraint*) constraint;
          if (intconstraint && 
              (lintval > intconstraint->max || lintval < intconstraint->min)) {
            manopt_help(optset, "unsigned int argument '%s' for option %c (%s) out of bounds\n",
              argset->args[arg].values[0], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
          }           
        }       
      }
      break;
    case INTRANGEOPT:
      if (noofvalues < 2) {
        manopt_help(optset, "range option %c (%s) requires at least two values",
            optset->opts[opt].shortopt, optset->opts[opt].longopt);
      }
      if (noofvalues > 2) {
        if(!lastarg) { 
          manopt_help(optset, "range option %c (%s) requires exactly two values",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
        } else { 
          manopt_unflag(argset, arg, 2);
        }
      } 
      if (!isint(argset->args[arg].values[0]) ||
          !isint(argset->args[arg].values[1]) ||
          (lintval=atoi(argset->args[arg].values[0])) == INT_MIN ||
          (rintval=atoi(argset->args[arg].values[1])) == INT_MIN ||
          lintval == INT_MAX || rintval == INT_MAX) {
        manopt_help(optset, "'%s'-'%s' for option %c (%s) out of range\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
 
      } else {
        if (lintval > rintval) {  
             manopt_help(optset, "'%s' > '%s' for option %c (%s)\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
 
        } else {
          intconstraint = (manopt_intconstraint*) constraint;
          if (intconstraint && 
              (rintval > intconstraint->max || lintval < intconstraint->min)) {
              manopt_help(optset, "'%s'-'%s' for option %c (%s) out of range\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
          } 
        }
      }       
      break;
    case UINTRANGEOPT:
      if (noofvalues < 2) {
        manopt_help(optset, "range option %c (%s) requires at least two values\n", 
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
      }
      if (noofvalues > 2) {
        if(!lastarg) { 
          manopt_help(optset,"range option %c (%s) requires exactly two values\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
        } else { 
          manopt_unflag(argset, arg, 2);
        }
      } 
      if (!isint(argset->args[arg].values[0]) ||
          !isint(argset->args[arg].values[1]) ||
          (lintval=atoi(argset->args[arg].values[0])) < 0 ||
          (rintval=atoi(argset->args[arg].values[1])) < 0 ||
          lintval == INT_MAX || rintval == INT_MAX) {
         manopt_help(optset, "'%s'-'%s' for option %c (%s) out of range\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
      } else {
        if (lintval > rintval) {  
              manopt_help(optset, "'%s'>'%s' for option %c (%s)\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
        } else {
          intconstraint = (manopt_intconstraint*) constraint;
          if (intconstraint && 
              (rintval > intconstraint->max || lintval < intconstraint->min)) {
            manopt_help(optset, "'%s'-'%s' for option %c (%s) out of range\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);
          } 
        }
      }       
      break;

    case DBLRANGEOPT:
      if (noofvalues < 2) {
         manopt_help(optset, "range option %c (%s) requires at least two values\n", 
             optset->opts[opt].shortopt, optset->opts[opt].longopt); 
 
      } else if (noofvalues > 2) {
        if(!lastarg) { 
           manopt_help(optset,"range option %c (%s) requires exactly two values\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt);  
        } else { 
          manopt_unflag(argset, arg, 2);
        }
      }
      if (!isfloat(argset->args[arg].values[0]) ||
          !isfloat(argset->args[arg].values[1]) ||
          (ldblval=atof(argset->args[arg].values[0])) == HUGE_VAL ||
          (rdblval=atof(argset->args[arg].values[1])) == HUGE_VAL) {
          manopt_help(optset, "'%s'-'%s' for option %c (%s) out of range\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);  
      } else {
        if (ldblval > rdblval) {
               manopt_help(optset, "'%s'>'%s' for option %c (%s)\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt);  
        } else {
          dblconstraint = (manopt_dblconstraint*) constraint;
          if (dblconstraint && 
              (rdblval > dblconstraint->max || ldblval < dblconstraint->min)) {
             manopt_help(optset, "'%s'-'%s' for option %c (%s) out of range\n",
              argset->args[arg].values[0], argset->args[arg].values[1], 
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
          } 
        }
      }      
      break; 
    case LISTOPT:
      if (noofvalues < 1) {
        manopt_help(optset, "list option %c (%s) requires at least one argument\n",
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
      } else {
        listconstraint = (manopt_listconstraint*) constraint;
        if (listconstraint) {
          if(noofvalues > listconstraint->maxlength) {
            if(!lastarg) {
              manopt_help(optset, "list option %c (%s) too long!",
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
            } else {
              manopt_unflag(argset, arg, listconstraint->maxlength);
            }
          }
          if(noofvalues < listconstraint->minlength) {
            if(!lastarg) {
              manopt_help(optset, "list option %c (%s) too short!",
              optset->opts[opt].shortopt, optset->opts[opt].longopt); 
            } else {
              manopt_unflag(argset, arg, listconstraint->maxlength);
            }
          }
        }
      }
      break;
    case SELECTOPT:
      if (noofvalues < 1) {
        manopt_help(optset, "list option %c (%s) requires at least one argument\n",
            optset->opts[opt].shortopt, optset->opts[opt].longopt);  
      } else {
        listconstraint = (manopt_listconstraint*) constraint;
        if(listconstraint) {
          if(noofvalues > listconstraint->maxlength) {
            if(!lastarg) {
              manopt_help(optset, "list option %c (%s) too long!",
                  optset->opts[opt].shortopt, optset->opts[opt].longopt);  
            } else { 
              manopt_unflag(argset, arg, listconstraint->maxlength);
            }
          } else if(noofvalues < listconstraint->minlength) {
            if(!lastarg) {
              manopt_help(optset, "list option %c (%s) too short!",
                  optset->opts[opt].shortopt, optset->opts[opt].longopt); 
            } else {
              manopt_unflag(argset, arg, listconstraint->maxlength);
            }
          } else {
            for(i=0; i < noofvalues; i++) {
              valid_select = (unsigned char) 0;
              for(j=0; j < listconstraint->noofitems; j++) { 
                if(!strcmp(listconstraint->items[j],
                      argset->args[arg].values[i])) {
                  valid_select = (unsigned char) 1;
                } 
              }
              if(!valid_select) {
                manopt_help(optset, "unknown value %s for select option %c (%s)",
                    argset->args[arg].values[i], 
                    optset->opts[opt].shortopt, optset->opts[opt].longopt);  

              }
            }
          }
        }
      }
      break;
    case MANOPT_BLOCKSEPARATOR:
      break;
    default:
      manopt_help(optset, "unkown option %s type\n", argset->args[arg].flagname);
      break;
  }

  return 1;
}


manopt_arg*
manopt_getopts(manopt_optionset* set, int argc, char **argv) {
  unsigned char *ucharval,
                optionfound = 0;
  unsigned int *uintval, *uintrangeval;
  char *charval, **ptr;
  int i, j, *intval, *intrangeval;
  double *dblval, *dblrangeval;
  manopt_argset argset;
  
  argset.noofargs = 0;
  argset.args = NULL;

  if(!manopt_parse_commandline(&argset, argc, argv)) {
    manopt_help(set, "error while parsing commandline.\n");
  }

  set->call = argv[0];

  for(j=0; j < argset.noofargs; j++) {
    if(argset.args[j].flagname) {
      optionfound = (unsigned char) 0;
      for(i=0; i < set->noofopts; i++) {
        if ((set->opts[i].longopt && 
              !strcmp(set->opts[i].longopt, argset.args[j].flagname)) 
            || (set->opts[i].shortopt && strlen(argset.args[j].flagname)==1 && 
              set->opts[i].shortopt == argset.args[j].flagname[0])) {
          if (set->opts[i].set) { 
            manopt_help(set, "option %c (%s) multiply selected!", 
                set->opts[i].longopt, set->opts[i].shortopt);
          } else { 
            set->opts[i].set = (unsigned char) 1;
            optionfound = 1;
            manopt_checkconstraint(set, i, &argset, j);
            memmove(&set->opts[i].arg, &argset.args[j], sizeof(manopt_arg));
            if(set->opts[i].reg_var) {
              switch(set->opts[i].type) {
                case FLAG:
                  ucharval = set->opts[i].reg_var;
                  *ucharval = 1;
                  break;
                case REQCHAROPT:
                case CHAROPT:
                  charval = set->opts[i].reg_var;
                  *charval = argset.args[j].values[0][0];
                  break;
                case REQUINTOPT:
                case REQINTOPT:
                  uintval = (unsigned int*) set->opts[i].reg_var;
                  *uintval = atoi(argset.args[j].values[0]);
                  break;
                case UINTOPT:
                case INTOPT:
                  intval = (int*) set->opts[i].reg_var;
                  *intval = atoi(argset.args[j].values[0]);
                  break;
                case REQDBLOPT:
                case DBLOPT:
                  dblval = (double*) set->opts[i].reg_var; 
                  *dblval = atof(argset.args[j].values[0]);
                  break;
                case REQSTRINGOPT:
                case STRINGOPT:
                  ptr = (char**) set->opts[i].reg_var;
                  ptr[0] = argset.args[j].values[0];
                  break;
                case INTRANGEOPT:
                  intrangeval = (int*) set->opts[i].reg_var; 
                  intrangeval[0] = atoi(argset.args[j].values[0]);
                  intrangeval[1] = atoi(argset.args[j].values[1]);
                  break;
                case UINTRANGEOPT:
                  uintrangeval = (unsigned int*) set->opts[i].reg_var;
                  uintrangeval[0] = atoi(argset.args[j].values[0]);
                  uintrangeval[1] = atoi(argset.args[j].values[1]);
                  break;
                case DBLRANGEOPT:
                  dblrangeval = (double*) set->opts[i].reg_var; 
                  dblrangeval[0] = atof(argset.args[j].values[0]);
                  dblrangeval[1] = atof(argset.args[j].values[1]);
                  break;
                default:
                  break;
              }
            }
          }
        }
      }
      if (!strcmp(argset.args[j].flagname,"h")||
          !strcmp(argset.args[j].flagname,"help")) {
        manopt_usage(set);
        exit(EXIT_FAILURE);
      }

      else if(!optionfound) {
        manopt_help(set, "option '%s' unknown\n", 
            argset.args[j].flagname); 
      }
    }
  } 
  argset.args = realloc(argset.args, sizeof(manopt_arg));

  for(i=0; i < set->noofopts; i++){
    if (set->opts[i].required && !set->opts[i].set) {
      manopt_help(set, "required option '%s' (%c) missing\n",
          set->opts[i].longopt, set->opts[i].shortopt);
    }
  }
  
  
  return &argset.args[0];
}

unsigned char
manopt_isset(manopt_optionset *set, char shortopt, char *longopt) {
  int i;
  
  for(i=0; i < set->noofopts; i++) {
    if((set->opts[i].shortopt == shortopt && set->opts[i].set) ||
        (set->opts[i].longopt && longopt &&
         !strcmp(set->opts[i].longopt,longopt) && set->opts[i].set )) {
      return 1;
    }
  }
  return 0;
}

manopt_arg*
manopt_getarg(manopt_optionset *set, char shortopt, char *longopt) {
  int i;
  
  for(i=0; i < set->noofopts; i++) {
    if((set->opts[i].shortopt == shortopt && set->opts[i].set) ||
        (set->opts[i].longopt && longopt &&
         !strcmp(set->opts[i].longopt,longopt) && set->opts[i].set )) {
      return &set->opts[i].arg;
    }
  }
  return NULL;
}



manopt_option*
manopt_longopt(manopt_optionset *set, char *longopt) {
  int i;
  
  for(i=0; i < set->noofopts; i++) {
    if(!strcmp(set->opts[i].longopt,longopt)) {
      return &set->opts[i];
    }
  }
  return NULL;
}

manopt_option*
manopt_shortopt(manopt_optionset *set, char shortopt){
int i;
  for(i=0; i < set->noofopts; i++) {
    if(set->opts[i].shortopt == shortopt) {
      return &set->opts[i];
    }
  }
  return NULL;
}


void
manopt_dumpoptionset(manopt_optionset *set) {
  int i,j;

  for(i=0; i < set->noofopts; i++) {
    printf("option: %s (%c)\n", set->opts[i].longopt, set->opts[i].shortopt);
    if(set->opts[i].arg.noofvalues) {
      for(j=0; j < set->opts[i].arg.noofvalues; j++) {
        printf("arg\n");
        printf("\t%s\n", set->opts[i].arg.values[j]);
      }
    }
  }
}

