/*
	  This file is part of gdub.
	  (C) 2006 Steve Hoffmann 
 
	  gdub is free software; you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation; either version 2, or (at your
	  option) any later version.
 
	  gdub is distributed in the hope that it will be useful, but
	  WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  General Public License for more details.
 
	  You should have received a copy of the GNU General Public License
	  along with gdub; see the file COPYING.  If not, write to the
	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
	  Boston, MA 02111-1307, USA.	
 
 */

/**
 * stack.c
 * implementation of a simple stack for int
 *
 * @author Steve Hoffmann
 * @email steve@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Tue Oct 28 10:42:34 CET 2008
 */

/*
 *  SVN
 *  Revision of last commit: $Rev: 73 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-10-29 10:03:28 +0100 (Wed, 29 Oct 2008) $
 *  Id: $Id: stack.c 73 2008-10-29 09:03:28Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/stack.c $
 */

#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include "stack.h"

/*----------------------------- bl_stackInit -----------------------------------
 *    
 * @brief 	init stack
 * @author 	Steve Hoffmann
 *   
 */
void bl_stackInit(Stack *stack, int allocelem) {
  if (allocelem <= 0){
    DBG("stack.c: Attempt to initialize a stack of size %d. Exit forced.\n",
	allocelem);
    exit(-1);
  }
  stack->stackspace = (Stackelement *) malloc(sizeof(Stackelement) * allocelem);
  if (stack->stackspace == NULL){
    DBG("stack.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }
  stack->allocelem=allocelem;
  stack->top=-1;

}

/*--------------------------- bl_stackDestruct ---------------------------------
 *    
 * @brief 	destruct stack
 * @author 	Steve Hoffmann
 *   
 */
void bl_stackDestruct(Stack *stack) {
  free(stack->stackspace);
  stack->top = 0;
  stack->allocelem = 0;
}

/*---------------------------- bl_stackIsEmpty ---------------------------------
 *    
 * @brief 	returns if the stack is empty
 * @author 	Steve Hoffmann
 *   
 */
BOOL bl_stackIsEmpty(Stack *stack) {
  return (stack->top < 0);
}


/*----------------------------- bl_stackPush -----------------------------------
 *    
 * @brief 	pushs elements on the top of the stack
 * @author 	Steve Hoffmann
 *   
 */		
void bl_stackPush(Stack *stack, Stackelement elem) {
		
  if(stack->top >= stack->allocelem - 1) {
			
    stack->stackspace = (Stackelement *) realloc(stack->stackspace, 
						 sizeof(Stackelement) * 
						 (stack->allocelem + BASEINC)); 
    if (stack->stackspace == NULL || BASEINC <= 0){
      DBG("stack.c: Memory reallocation failed. Exit forced.\n", NULL);
      exit(-1);
    }
    stack->allocelem += BASEINC; 
  }

  stack->top++;
  stack->stackspace[stack->top] = elem;
}

/*------------------------------ bl_stackPop -----------------------------------
 *    
 * @brief 	pops the top of the stack
 * @author 	Steve Hoffmann
 *   
 */
Stackelement bl_stackPop(Stack *stack){
  if(!bl_stackIsEmpty(stack)) {
    return stack->stackspace[stack->top--];	
  }	
  return STACK_NULL_TYPE;
}

/*------------------------------ bl_stackTop -----------------------------------
 *    
 * @brief 	returns top of the stack
 * @author 	Steve Hoffmann
 *   
 */
Stackelement bl_stackTop(Stack *stack){	
  if(!bl_stackIsEmpty(stack)) {
    return stack->stackspace[stack->top];
  }
  return STACK_NULL_TYPE;
}
 
/*------------------------------ bl_stackTopN --------------------------------
 *    
 * @brief 	returns Nth highest object of the stack
 *              with N = 0,..,numofelems - 1
 * @author 	Steve Hoffmann
 *   
 */
Stackelement bl_stackTopN(Stack *stack, Uint n){	
  if(!bl_stackIsEmpty(stack) && n >= 0 && n <= stack->top) {
    return stack->stackspace[stack->top - n];
  }
  return STACK_NULL_TYPE;
}

/*------------------------------ bl_stackSize ----------------------------------
 *    
 * @brief 	returns number of elements on the stack
 * @author 	Steve Hoffmann
 *   
 */
Uint bl_stackSize(Stack *stack) {
  return (stack->top + 1);    
}

