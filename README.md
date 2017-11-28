# barebones_C_miniapp_CoreNeuron_kernels

This miniapp was developed by Francesco Cremonesi as a barebones version of the more complete 
[BlueBrain/neuromapp](https://github.com/BlueBrain/neuromapp.git) application.

The purpose of the original Neuronm(ini)app(lication)library is to reproduce the algorithms of the 
main software of the BBP as a collection of mini-apps.

The barebones_C_miniapp_CoreNeuron_kernels instead is solely focused on the synapse kernels.
The purpose of this duplication is to avoid any overhead coming form the original BlueBrain/neuromapp
infrastructure, and to allow the author to change and instrument the kernel code freely.

## BUILD

The barebones_C_miniapp_CoreNeuron_kernels depends on `openmp` for shared-memory parallelization. 
Make sure you can compile openmp code and, if necessary, edit the Makefile to provide the correct flags to the compiler.

To build barebones_C_miniapp_CoreNeuron_kernels, edit the Makefile with the correct compiler and flags combination, 
then type `make` in your command line.

#### Using the likwid toolsuite

The barebones_C_miniapp_CoreNeuron_kernels is already instrumented to allow for benchmarking and profiling using the
[likwid marker API](https://github.com/RRZE-HPC/likwid/wiki/TutorialMarkerC).
To enable the likwid marker API, build the `likwid` target instead by typing `make likwid` in your command line.

#### Selecting different kernels

By default, the `current` function of the `ProbAMPANMDA_EMS` synapse is selected. If you want to change this behaviour, 
comment/uncomment the relevant lines in the `main.c` file:
```c
    for(i=0 ; i < NDUP; ++i) {
        for(j=0 ; j < NITER; ++j) {
//            state_Exp2Syn(ntu[i] ,mech_id);
//            current_Exp2Syn(ntu[i], mech_id);
//            state_ProbAMPANMDA_EMS(ntu[i], mech_id);
             current_ProbAMPANMDA_EMS(ntu[i], mech_id);
}
}
```
