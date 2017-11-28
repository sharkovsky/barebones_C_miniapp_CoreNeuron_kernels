/* Created by Language version: 6.2.0 */
/* VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <omp.h>

#include "common/memory/nrnthread.h"
#include "exp2syn.h"

#undef PI
#define nil 0

#define STRIDE _cntml + _iml

#define t _nt->_t
#define dt _nt->_dt
#define tau1 _p[0*STRIDE]
#define tau2 _p[1*STRIDE]
#define e _p[2*STRIDE]
#define i _p[3*STRIDE]
#define g _p[4*STRIDE]
#define A _p[5*STRIDE]
#define B _p[6*STRIDE]
#define factor _p[7*STRIDE]
#define DA _p[8*STRIDE]
#define DB _p[9*STRIDE]
#define _v_unused _p[10*STRIDE]
#define _g_unused _p[11*STRIDE]
#define _tsav _p[12*STRIDE]
#define _nd_area  _nt_data[_ppvar[0*STRIDE]]

#define initmodel initmodel_Exp2Syn
#define current current_Exp2Syn
#define state state_Exp2Syn
#define output_states output_states_Exp2Syn


void initmodel(NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml;

    // #pragma omp for
    for (_iml = 0; _iml < _cntml; ++_iml) {
        if ( tau1 / tau2 > 0.9999 ) {
            tau1 = 0.9999 * tau2 ;
        }
        if ( tau1 / tau2 < 1e-9 ) {
            tau1 = tau2 * 1e-9 ;
        }
        A = 0.0 ;
        B = 0.0 ;
    }
}


void current(NrnThread* _nt, int mech_id){
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;

    int _iml, _nd_idx;
    double _v, _lrhs, _lg, _mfact;

    double * restrict _p = _ml.data;
    // TODO: check if I can actually restrict here, and ask why the origi    nal nmodl vectorizes
    double * restrict _vec_shadow_rhs = _nt->_shadow_rhs;
    double * restrict _vec_shadow_d = _nt->_shadow_d;
    double * _vec_rhs = _nt->_actual_rhs;
    double * _vec_d = _nt->_actual_d;
    double * _vec_v = _nt->_actual_v;
    double * _nt_data = _nt->_data;

    int * _ppvar;
    int * _ni = _ml.nodeindices;
    _ppvar = _ml.pdata;

    for (_iml = 0; _iml < _cntml; ++_iml) {
        _nd_idx = _ni[_iml];
        _v = _vec_v[_nd_idx];

        g = B - A ;
        i = g * ( _v - e ) ;

        _lrhs = i;
        _lg = g;
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

int  state (NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml;

    // #pragma omp for
    for (_iml = 0; _iml < _cntml; ++_iml) {
        A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
        B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
    }
    return 0;
}

void output_states(NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml;


    for (_iml = 0; _iml < _cntml; ++_iml) {
        printf("A: %g\n", A);
        printf("B: %g\n", B);
    }
}

