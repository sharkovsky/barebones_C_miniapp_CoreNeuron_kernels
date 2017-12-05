# barebones_C_miniapp_CoreNeuron_kernels

This miniapp was developed by Francesco Cremonesi as a barebones version of the more complete 
[BlueBrain/neuromapp](https://github.com/BlueBrain/neuromapp.git) application.

The purpose of the original Neuronm(ini)app(lication)library is to reproduce the algorithms of the 
main software of the BBP as a collection of mini-apps.

The barebones_C_miniapp_CoreNeuron_kernels instead is solely focused on the synapse kernels.
The purpose of this duplication is to avoid any overhead coming form the original BlueBrain/neuromapp
infrastructure, and to allow the author to change and instrument the kernel code freely.

## BUILD

To build barebones_C_miniapp_CoreNeuron_kernels, edit the Makefile with the correct compiler and flags combination, 
then type `make` in your command line.

#### openmp dependency

The barebones_C_miniapp_CoreNeuron_kernels depends on `openmp` for shared-memory parallelization. 
Make sure you can compile openmp code and, if necessary, edit the Makefile to provide the correct flags to the compiler.

#### Using the likwid toolsuite

The barebones_C_miniapp_CoreNeuron_kernels is already instrumented to allow for benchmarking and profiling using the
[likwid marker API](https://github.com/RRZE-HPC/likwid/wiki/TutorialMarkerC) from the [RRZE-HPC/likwid](https://github.com/RRZE-HPC/likwid) toolsuite.
To enable the likwid marker API, build the `likwid` target instead by typing `make likwid` in your command line.

#### Selecting different kernels

For now, two synapse types are supported:
1. Exp2Syn is a simple static synapse model with only two state variables;
2. ProbAMPANMDA_EMS is a more complex synapse model with short-term plasticity dynamics, that is modeled using four state variables.

By default, the `current` function of the `ProbAMPANMDA_EMS` synapse is selected. If you want to change this behaviour, you must do two things:
1. comment/uncomment the relevant lines in the `main.c` file
```c
    for(i=0 ; i < NDUP; ++i) {
        for(j=0 ; j < NITER; ++j) {
#ifdef PROBAMPA
//        state_ProbAMPANMDA_EMS(ntu[i], mech_id);
        current_ProbAMPANMDA_EMS(ntu[i], mech_id);
#elif defined EXP2SYN
//        state_Exp2Syn(ntu[i] ,mech_id);
//        current_Exp2Syn(ntu[i], mech_id);
#endif
        }
    }
```
2. rerun make and, if necessary, define the correct value for `SYN_FLAG` as such:
```bash
make SYN_FLAG=-DEXP2SYN
```

## RUN

The barebones_C_miniapp_CoreNeuron_kernels takes three arguments as input:

- the input data file ( usually named `*.nrnthread.txt`);
- the number of data duplications to be stored internally;
- the number of _time_ iterations, or iterations of the kernel loop;

```bash
./kernels_app.exe <input file> <n times data is duplicated> <n iter main loop>
e.g. ./kernels_app.exe nrnthread.txt 1 100
```

#### example dataset

Datasets for the miniapp are difficult to generate.
Among the datasets that come with this distribution, a very small example one is given.
You can test that the BUILD process executed correctly by running and validating the output given below
```bash
$ ./kernels_app.exe data/example.probampa.nrnthread.txts 1 1
num_omp_threads 1
Reading pdata because szdp = 3
Finished reading file
Found index of Mechanism type 53 at id 1
Number of instances: 128
Finished initialization
Executing 1 iterations on 1 data duplications
Finished computing
```

#### likwid benchmarking and profiling

The barebones_C_miniapp_CoreNeuron_kernels already has the necessary code instrumentation to support profiling using the likwid toolsuite.
To profile a given kernel, make sure that you have compiled following the instructions for likwid given in the BUILD section.
Then, you can run it using the `likwid-perfctr` tool as follows:
```bash
likwid-perfctr -C <your core bindings> -g <your desired analysis> -m -s 0x0
```
For more information on how to use `likwid-perfctr`, please consult the [wiki](https://github.com/RRZE-HPC/likwid/wiki/likwid-perfctr) page.
