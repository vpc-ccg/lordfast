CC          ?= gcc
WRAP_MALLOC  = -DUSE_MALLOC_WRAPPERS
DFLAGS       = -DHAVE_PTHREAD $(WRAP_MALLOC)
CFLAGS       = -Wall -Wno-unused-function -Wno-unused-but-set-variable -O2 $(DFLAGS)
BWAOBJ       = utils.o malloc_wrap.o bntseq.o bwt_gen.o bwt.o bwtindex.o bwa.o is.o QSufSort.o ksw.o

all: $(BWAOBJ)

clean:
	@rm -f $(BWAOBJ)
