## GCC
#CC=gcc
#CFLAGS=-O2 -Wno-unused-result -g -ftree-vectorize -fopt-info-vec-optimized -DLIKWID_PERFMON -fopenmp
#INCFLAGS=-I. $(LIKWID_INC)
#LINKFLAGS=-lm -fopenmp $(LIKWID_LIB) -llikwid

## INTEL
CC=icc
CFLAGS=-O2 -g -qopt-report-phase=vec -qopt-streaming-stores never -restrict -unroll0 -qopenmp -DLIKWID_PERFMON
INCFLAGS=-I. $(LIKWID_INC)
LINKFLAGS=-lm -qopenmp $(LIKWID_LIB) -llikwid

.PHONY: all clean

all: kernels_app.exe

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

