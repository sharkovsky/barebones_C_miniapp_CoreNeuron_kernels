SYN_FLAG=-DPROBAMPA
#SYN_FLAG=-DEXP2SYN

## GCC
#CC=gcc
#CFLAGS=-O2 -Wno-unused-result -g -ftree-vectorize -fopt-info-vec-optimized -DLIKWID_PERFMON -fopenmp $(SYN_FLAG)
#INCFLAGS=-I. $(LIKWID_INC)
#LINKFLAGS=-lm -fopenmp $(LIKWID_LIB) -llikwid

## INTEL
CC=icc
CFLAGS=-O2 -g -qopt-report-phase=vec -qopt-streaming-stores never -restrict -unroll0 -qopenmp $(SYN_FLAG)
INCFLAGS=-I.
LINKFLAGS=-lm -qopenmp

.PHONY: all clean likwid

all: kernels_app.exe

likwid: CFLAGS+=-DLIKWID_PERFMON

likwid:INCFLAGS+=$(LIKWID_INC)

likwid:LINKFLAGS+=$(LIKWID_LIB) -llikwid

likwid: all

kernels_app.exe: main.o ProbAMPANMDA_EMS.o exp2syn.o nrnthread.o memory.o
	$(CC) -o $@ $^ $(LINKFLAGS)

ProbAMPANMDA_EMS.o: ProbAMPANMDA_EMS.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $^

exp2syn.o: exp2syn.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $^

main.o: main.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $^

nrnthread.o: common/memory/nrnthread.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $^

memory.o: common/memory/memory.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $^

clean:
	rm -f *.o *.optrpt *.s *.pp

