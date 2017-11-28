/*
 * Neuromapp - main.c, Copyright (c), 2015,
 * Timothee Ewart - Swiss Federal Institute of technology in Lausanne,
 * Pramod Kumbhar - Swiss Federal Institute of technology in Lausanne,
 * timothee.ewart@epfl.ch,
 * paramod.kumbhar@epfl.ch
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.
 */

/**
 * @file neuromapp/coreneuron_1.0/cstep/main.c
 * \brief Implements a miniapp combining  the compute kernel and the Hines solver of coreneuron 1.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <omp.h>

#include "common/memory/nrnthread.h"

#include "ProbAMPANMDA_EMS.h"
#include "exp2syn.h"

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

// Get OMP header if available

/** \fn compute_wrapper(NrnThread *nt, struct input_parameters* p)
    \brief Start the computation of kernel following the input parameter
    \param nt the data structure where all the datas are saved
    \param p input parameters where are defined the wanted computation
 */

void compute_wrapper(NrnThread *nt, int mech_id, char syn_type);

int main(int argc, char *const argv[]) {

    if (argc != 5) {
        printf("Wrong number of arguments! Usage:\n");
        printf("%s <input file> <n times data is duplicated> <n iter main loop> <synapse type: p for ProbAMPA or e for Exp2Syn >\n", argv[0]);
        printf("e.g. %s nrnthread.txt 1 100 a\n", argv[0]);
        return -1;
    }


    int NDUP = atoi(argv[2]);
    int NITER = atoi(argv[3]);
    int i,j;

    LIKWID_MARKER_INIT;

#pragma omp parallel
    {
        LIKWID_MARKER_THREADINIT;
    }

    #pragma omp parallel
    {
        if ( omp_get_thread_num() == 0 ) {
        printf("num_omp_threads %d\n", omp_get_num_threads());
        }
    }


    int mech_type_to_find;

    char c = argv[4][0];

    // select the synapse type
    switch( c ) {
        case 'p' :
            mech_type_to_find = 53;
            break;
        case 'e' :
            mech_type_to_find = 10;
            break;
        default :
            printf("error: unsupported synapses type %s, supported types are p,e\n", argv[4]);
            return 1;
            break;
    }

    //duplicate the data
    NrnThread * ntu[NDUP];
    FILE * fp;
    int flag_file_problem = 0;

#pragma omp parallel
    {
        LIKWID_MARKER_START("readInput");
#pragma omp for private(fp) schedule(static, 1)
    for(i=0; i<NDUP; ++i ) {

        fp = fopen( argv[1], "r");
        if ( fp == NULL ) {
            printf("error opening file %s\n", argv[1]);
            flag_file_problem = 1;
        }

        ntu[i] = (NrnThread *) malloc( sizeof(NrnThread) );
        //nrnthread_copy(nt, ntu[i]);
        nrnthread_read(fp, ntu[i]);
        fclose(fp);
    }
        LIKWID_MARKER_STOP("readInput");
    }
    printf("Finished reading file\n");

    if ( flag_file_problem ) {
        return 1;
    }

    // find index of the synapse
   int mech_id;
   int flag = 0;
    for(mech_id=0 ; mech_id < ntu[0]->nmech; ++mech_id){
        if ( ntu[0]->ml[mech_id].type == mech_type_to_find ) {
            printf("Found index of Mechanism type %d at id %d\n", mech_type_to_find, mech_id);
            printf("Number of instances: %d\n", ntu[0]->ml[mech_id].nodecount);
            flag = 1;
            break;
        }
    }

    if ( flag == 0 ){
        printf("Error: could not find mechanism type %d. Exiting\n", mech_type_to_find);
        return 2;
    }

#pragma omp parallel for schedule(static, 1)
    for(i=0; i<NDUP; ++i ) {
        switch( c ) {
            case 'p' :
                initmodel_ProbAMPANMDA_EMS(ntu[i], mech_id);
                break;
            case 'e' :
                initmodel_Exp2Syn(ntu[i], mech_id);
                break;
        }
    }

    //output_states(ntu[0], mech_id);
    printf("Finished initialization\n");
    printf("Executing %d iterations on %d data duplications\n", NITER, NDUP);
    fflush(stdout);

#pragma omp parallel
    {
        LIKWID_MARKER_START("kernelloop");
#pragma omp for private(j) schedule(static, 1)
    for(i=0 ; i < NDUP; ++i) {
        for(j=0 ; j < NITER; ++j) {
//            state_Exp2Syn(ntu[i] ,mech_id);
//            current_Exp2Syn(ntu[i], mech_id);
//            state_ProbAMPANMDA_EMS(ntu[i], mech_id);
            current_ProbAMPANMDA_EMS(ntu[i], mech_id);
        }
    }
        LIKWID_MARKER_STOP("kernelloop");
    }

    printf("Finished computing\n");
    fflush(stdout);
    //output_states(ntu[0], mech_id);
    //fflush(stdout);

    for(i=0 ; i < NDUP; ++i)
        nrnthread_dealloc(ntu[i]);

    LIKWID_MARKER_CLOSE;
    return 0;
}













void compute_wrapper(NrnThread *nt, int mech_id, char syn_type) {
//    switch( syn_type ) {
//        case 'p' :
//            state_ProbAMPANMDA_EMS(nt, mech_id);
//            current_ProbAMPANMDA_EMS(nt, mech_id);
//            break;
//        case 'e' :
//            state_Exp2Syn(nt ,mech_id);
//            current_Exp2Syn(nt, mech_id);
//            break;
//    }
}
