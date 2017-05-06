CC 			?= gcc

INCS		:= 
LIBS		:= -lm -lssl -lpthread
CFLAGS 		:= -w -pedantic -std=c99 -O3 -DPROG3NFO -DDBGLEVEL=0
LDFLAGS 	:= 
# CFLAGS=-Wall -pedantic -std=c99 -O3 -DPROG3NFO -DDBGLEVEL=0
# LDFLAGS=-lm -lssl -lpthread
# LIBS=-lob -lm -lopenssl -lpthread

CLASPOBJ 	 = stringutils.o \
			   mathematics.o \
			   memman.o \
			   fileio.o \
			   debug.o \
			   info.o \
			   stack.o \
			   sort.o \
			   list.o \
			   vtprogressbar.o \
			   queue.o \
			   manopt.o \
			   sltypes.o \
			   slchain.o \
			   bintree.o \
			   rangetree.o \
			   vebtree.o \
			   container.o \
			   vstack.o \
			   vqueue.o \

.PHONY: clasp

clasp: ${CLASPOBJ}

clean: 
	@rm -f ${CLASPOBJ}
