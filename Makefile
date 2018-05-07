PROG_VERSION := 0.0.6
VERBOSITY := 0

all: OPTIMIZE_FLAGS build
log1: LOG1_FLAG OPTIMIZE_FLAGS build
log2: LOG2_FLAG OPTIMIZE_FLAGS build
log3: LOG3_FLAG OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build: clean-exe SSE_FLAGS clasplib bwalib lordfast clean

CC          ?= gcc
CXX         ?= g++

SRCDIR      := src
CLASPDIR    := $(SRCDIR)/clasp
BWADIR      := $(SRCDIR)/bwa
EDLIBDIR    := $(SRCDIR)/edlib

INCS        := -I$(CLASPDIR) -I$(BWADIR) -I$(EDLIBDIR)
LIBS        := -lz -lm -lpthread
CFLAGS      := -w -fno-pic
CXXFLAGS    := -w -fno-pic -DPROG_VERSION=\"$(PROG_VERSION)\" $(INCS)
LDFLAGS     := #-static

LORDFASTOBJ   = $(SRCDIR)/baseFAST.o \
               $(SRCDIR)/Common.o \
               $(SRCDIR)/CommandLineParser.o \
               $(SRCDIR)/Reads.o \
               $(SRCDIR)/BWT.o \
               $(SRCDIR)/LordFAST.o \
               $(SRCDIR)/Chain.o \
               $(EDLIBDIR)/edlib.o \
               HELP.o \

CLASPOBJ     = $(CLASPDIR)/*.o

BWAOBJ       = $(BWADIR)/*.o

lordfast: $(LORDFASTOBJ)
ifeq ($(shell uname -s),Linux)
	$(CXX) -w $(LORDFASTOBJ) $(BWAOBJ) $(CLASPOBJ) -o $@ ${LDFLAGS} ${LIBS}
else
	$(CXX) -Wl,-no_pie -fno-pic -w $(LORDFASTOBJ) $(BWAOBJ) $(CLASPOBJ) -o $@ ${LDFLAGS} ${LIBS}
endif

clasplib: 
	@$(MAKE) -C $(CLASPDIR)

bwalib:
	@$(MAKE) -C $(BWADIR)

clean:
	@rm -f $(LORDFASTOBJ)
	@rm -f HELPstub.c
	@rm -f HELPstub.o

clean-exe:
	@rm -f lordfast
	@rm -f HELP

clean-all: clean-exe clean
	@$(MAKE) -C $(CLASPDIR) clean
	@$(MAKE) -C $(BWADIR) clean

HELP.o:
	@groff -Tascii -man HELP.man > HELP
ifeq ($(shell uname -s),Linux)
	@ld -r -b binary -o HELP.o HELP
else
	@touch HELPstub.c
	gcc -o HELPstub.o -c HELPstub.c
	ld -r -o HELP.o -sectcreate binary HELP HELP HELPstub.o
endif

DEBUG_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -O2 -DVERBOSITY=$(VERBOSITY))

LOG1_FLAG:
	$(eval VERBOSITY = 1)

LOG2_FLAG:
	$(eval VERBOSITY = 2)

LOG3_FLAG:
	$(eval VERBOSITY = 3)

PROFILE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)

SSE_FLAGS:
ifeq ($(shell uname -s),Linux)
ifeq ($(with-sse4),no)
		$(shell echo "-DSSE4=0")
else
        	$(eval CXXFLAGS = $(CXXFLAGS) \
        	$(shell gv=`gcc -dumpversion`; \
            	    sc=`grep -c "sse4" /proc/cpuinfo`; \
                	echo $$sc.$$gv | awk -F. '{if($$1>0 && ($$2>=4 || ($$2>=4 && $$3>=4))) print "-DSSE4=1 -msse4.2"; else print "-DSSE4=0"}'))
endif
else
ifeq ($(with-sse4),no)
		$(shell echo "-DSSE4=0")
else
        $(eval CXXFLAGS = $(CXXFLAGS) \
        $(shell gv=`gcc -dumpversion`; \
                sc=`sysctl -n machdep.cpu.features | grep -c "SSE4"` ;\
                echo $$sc.$$gv | awk -F. '{if($$1>0 && ($$2>=4 || ($$2>=4 && $$3>=4))) print "-DSSE4=1 -msse4.2"; else print "-DSSE4=0"}'))
endif
endif
