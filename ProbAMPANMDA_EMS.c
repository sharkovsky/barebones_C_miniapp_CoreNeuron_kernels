/* Created by Language version: 7.5.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include "myexp.h"

#include "common/memory/nrnthread.h"
#include "ProbAMPANMDA_EMS.h"

#undef PI
#define nil 0

#define STRIDE _cntml + _iml

#define t _nt->_t
#define dt _nt->_dt
#define tau_r_AMPA _p[0*STRIDE]
#define tau_d_AMPA _p[1*STRIDE]
#define tau_r_NMDA _p[2*STRIDE]
#define tau_d_NMDA _p[3*STRIDE]
#define Use _p[4*STRIDE]
#define Dep _p[5*STRIDE]
#define Fac _p[6*STRIDE]
#define e _p[7*STRIDE]
#define mg _p[8*STRIDE]
#define u0 _p[9*STRIDE]
#define Nrrp _p[10*STRIDE]
#define synapseID _p[11*STRIDE]
#define verboseLevel _p[12*STRIDE]
#define NMDA_ratio _p[13*STRIDE]
#define i _p[14*STRIDE]
#define i_AMPA _p[15*STRIDE]
#define i_NMDA _p[16*STRIDE]
#define g_AMPA _p[17*STRIDE]
#define g_NMDA _p[18*STRIDE]
#define g _p[19*STRIDE]
#define unoccupied _p[20*STRIDE]
#define occupied _p[21*STRIDE]
#define tsyn _p[22*STRIDE]
#define u _p[23*STRIDE]
#define A_AMPA _p[24*STRIDE]
#define B_AMPA _p[25*STRIDE]
#define A_NMDA _p[26*STRIDE]
#define B_NMDA _p[27*STRIDE]
#define factor_AMPA _p[28*STRIDE]
#define factor_NMDA _p[29*STRIDE]
#define mggate _p[30*STRIDE]
#define DA_AMPA _p[31*STRIDE]
#define DB_AMPA _p[32*STRIDE]
#define DA_NMDA _p[33*STRIDE]
#define DB_NMDA _p[34*STRIDE]
#define v _p[35*STRIDE]
#define _g _p[36*STRIDE]
#define _tsav _p[37*STRIDE]
#define _nd_area  _nt_data[_ppvar[0*STRIDE]]

#define initmodel initmodel_ProbAMPANMDA_EMS
#define current current_ProbAMPANMDA_EMS
#define state state_ProbAMPANMDA_EMS
#define output_states output_states_ProbAMPANMDA_EMS

void initmodel(NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml;

    // #pragma omp for
    for (_iml = 0; _iml < _cntml; ++_iml) {
        tau_r_AMPA = 0.2;
        tau_d_AMPA = 1.7;
        tau_r_NMDA = 0.29;
        tau_d_NMDA = 43;
        Use = 1;
        Dep = 100;
        Fac = 10;
        e = 0;
        mg = 1;
        u0 = 0;
        Nrrp = 1;
        synapseID = 0;
        verboseLevel = 0;
        NMDA_ratio = 0.71;
        A_NMDA = 10.;
        A_AMPA = 10.;
        B_NMDA = 10.;
        B_AMPA = 10.;
        tsyn = 0.0 ;
        u = 1.0 ;
        unoccupied = 0.0 ;
        occupied = Nrrp ;
    }
}

int  state (NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml;

    #pragma omp simd
    for (_iml = 0; _iml < _cntml; ++_iml) {
    A_AMPA = A_AMPA + (1. - exp(dt*(( 1.0 ) / tau_r_AMPA)))*(- ( 0.0 ) / ( ( 1.0 ) / tau_r_AMPA ) - A_AMPA) ;
    B_AMPA = B_AMPA + (1. - exp(dt*(( 1.0 ) / tau_d_AMPA)))*(- ( 0.0 ) / ( ( 1.0 ) / tau_d_AMPA ) - B_AMPA) ;
    A_NMDA = A_NMDA + (1. - exp(dt*(( 1.0 ) / tau_r_NMDA)))*(- ( 0.0 ) / ( ( 1.0 ) / tau_r_NMDA ) - A_NMDA) ;
    B_NMDA = B_NMDA + (1. - exp(dt*(( 1.0 ) / tau_d_NMDA)))*(- ( 0.0 ) / ( ( 1.0 ) / tau_d_NMDA ) - B_NMDA) ;
    }
    return 0;
}

void current(NrnThread* _nt, int mech_id){
    double gmax = 0.001;

    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;

    int _iml, _nd_idx;
    double _v, _lrhs, _lg, _mfact;

    double * restrict _p = _ml.data;
    // TODO: check if I can actually restrict here, and ask why the original nmodl vectorizes
    double * restrict _vec_shadow_rhs = _nt->_shadow_rhs;
    double * restrict _vec_shadow_d = _nt->_shadow_d;
    double * _vec_rhs = _nt->_actual_rhs;
    double * _vec_d = _nt->_actual_d;
    double * _vec_v = _nt->_actual_v;
    double * _nt_data = _nt->_data;

    int * _ppvar;
    int * _ni = _ml.nodeindices;
    _ppvar = _ml.pdata;

#pragma omp simd
    for (_iml = 0; _iml < _cntml; ++_iml) {
        _nd_idx = _ni[_iml];
        _v = _vec_v[_nd_idx];

        mggate = 1.0 / ( 1.0 + exp ( 0.062 * - ( _v ) ) * ( mg / 3.57 ) ) ;
        g_AMPA = gmax * ( B_AMPA - A_AMPA ) ;
        g_NMDA = gmax * ( B_NMDA - A_NMDA ) * mggate ;
        g = g_AMPA + g_NMDA ;
        i_AMPA = g_AMPA * ( _v - e ) ;
        i_NMDA = g_NMDA * ( _v - e ) ;
        i = i_AMPA + i_NMDA ;

        _lrhs = i;
        _lg = g_AMPA + g_NMDA;
        _mfact =  1.e2/(_nd_area);
        _lg *=  _mfact;
        _lrhs *= _mfact;

        _vec_shadow_rhs[_iml] = _lrhs;
        _vec_shadow_d[_iml] = _lg;

    }

    /*
    for (_iml = 0; _iml < _cntml; ++_iml) {
      _nd_idx = _ni[_iml];
      _vec_rhs[_nd_idx] -= _vec_shadow_rhs[_iml];
      _vec_d[_nd_idx] += _vec_shadow_d[_iml];
    }
    */

}

void output_states(NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml;

    printf("nodecount: %d\n", _cntml);

    for (_iml = 0; _iml < _cntml; ++_iml) {
        printf("A_AMPA: %g\n", A_AMPA);
        printf("B_AMPA: %g\n", B_AMPA);
        printf("A_NMDA: %g\n", A_NMDA);
        printf("B_NMDA: %g\n", B_NMDA);
    }
}

