/* Created by Language version: 7.5.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include "common/memory/nrnthread.h"
#include "NaTs2_t.h"

#undef PI
#define nil 0

#define _STRIDE _cntml + _iml

#define t _nt->_t
#define dt _nt->_dt
#define gNaTs2_tbar _p[0*_STRIDE]
#define ina _p[1*_STRIDE]
#define gNaTs2_t _p[2*_STRIDE]
#define m _p[3*_STRIDE]
#define h _p[4*_STRIDE]
#define ena _p[5*_STRIDE]
#define mInf _p[6*_STRIDE]
#define mTau _p[7*_STRIDE]
#define mAlpha _p[8*_STRIDE]
#define mBeta _p[9*_STRIDE]
#define hInf _p[10*_STRIDE]
#define hTau _p[11*_STRIDE]
#define hAlpha _p[12*_STRIDE]
#define hBeta _p[13*_STRIDE]
#define Dm _p[14*_STRIDE]
#define Dh _p[15*_STRIDE]
#define _v_unused _p[16*_STRIDE]
#define _g_unused _p[17*_STRIDE]

#define _ion_ena	_nt_data[_ppvar[0*_STRIDE]]
#define _ion_ina	_nt_data[_ppvar[1*_STRIDE]]
#define _ion_dinadv	_nt_data[_ppvar[2*_STRIDE]]

#define current current_NaTs2_t
#define state state_NaTs2_t
#define initmodel initmodel_NaTs2_t
#define output_states output_states_NaTs2_t

int  state (NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml;
    double _lqt ;

    int * _ni = _ml.nodeindices;
    double * _vec_v = _nt->_actual_v;

    int _nd_idx;
    double v;

    #pragma omp simd
    for (_iml = 0; _iml < _cntml; ++_iml) {
        _nd_idx = _ni[_iml];
        v = _vec_v[_nd_idx];
		mAlpha = ( 0.182 * ( v + 32.0 ) ) / ( 1.0 - ( exp ( - ( v + 32.0 ) / 6.0 ) ) ) ;
		mBeta = ( 0.124 * ( - v - 32.0 ) ) / ( 1.0 - ( exp ( - ( - v - 32.0 ) / 6.0 ) ) );
		mInf = mAlpha / ( mAlpha + mBeta ) ;
		mTau = ( 1.0 / ( mAlpha + mBeta ) ) / 2.95 ;
		hAlpha = ( - 0.015 * ( v + 60.0 ) ) / ( 1.0 - ( exp ( ( v + 60.0 ) / 6.0 ) ) ) ;
		hBeta = ( - 0.015 * ( - v - 60.0 ) ) / ( 1.0 - ( exp ( ( - v - 60.0 ) / 6.0 ) ) );
		hInf = hAlpha / ( hAlpha + hBeta ) ;
		hTau = ( 1.0 / ( hAlpha + hBeta ) ) / 2.95 ;
		m = m + (1. - exp(dt*( - 1.0 / mTau)))*(- ( mInf / mTau ) / ( - 1.0 / mTau ) - m) ;
		h = h + (1. - exp(dt*( - 1.0 / hTau)))*(- ( hInf / hTau ) / ( - 1.0 / hTau ) - h) ;
	}
return 0;
}

void current(NrnThread* _nt, int mech_id){
    double gmax = 0.001;

    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;

    int _iml, _nd_idx;
    double _v, _rhs, _g, _mfact;

    double * restrict _p = _ml.data;
    // TODO: check if I can actually restrict here, and ask why the original nmodl vectorizes
    double * _vec_rhs = _nt->_actual_rhs;
    double * _vec_d = _nt->_actual_d;
    double * _vec_v = _nt->_actual_v;
    double * _nt_data = _nt->_data;

    int * _ppvar;
    int * _ni = _ml.nodeindices;
    _ppvar = _ml.pdata;

#pragma ivdep
#pragma vector aligned
#pragma vector always
    for (_iml = 0; _iml < _cntml; ++_iml) {
        _nd_idx = _ni[_iml];
        _v = _vec_v[_nd_idx];
        ena = _ion_ena;
        gNaTs2_t = gNaTs2_tbar * m * m * m * h ;
        ina = gNaTs2_t * ( _v - ena ) ;
        _rhs = ina;
        _g = gNaTs2_t;
        _ion_ina += ina ;
#pragma distribute_point
        _vec_rhs[_nd_idx] -= _rhs;
        _vec_d[_nd_idx] += _g;

    }

}

void initmodel(NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml, _nd_idx;
    double _v;
    int * _ni = _ml.nodeindices;
    int * _ppvar;
    _ppvar = _ml.pdata;
    double * _vec_v = _nt->_actual_v;
    double * _nt_data = _nt->_data;

    // #pragma omp for
    for (_iml = 0; _iml < _cntml; ++_iml) {
        _nd_idx = _ni[_iml];
        _v = _vec_v[_nd_idx];
        gNaTs2_tbar = 0.00001;
        _ion_ena = 60.;
        _ion_ina = 0.;

		mAlpha = ( 0.182 * ( _v + 32.0 ) ) / ( 1.0 - ( exp ( - ( _v + 32.0 ) / 6.0 ) ) ) ;
		mBeta = ( 0.124 * ( - _v - 32.0 ) ) / ( 1.0 - ( exp ( - ( - _v - 32.0 ) / 6.0 ) ) );
		mInf = mAlpha / ( mAlpha + mBeta ) ;
        m = mInf;

		hAlpha = ( - 0.015 * ( _v + 60.0 ) ) / ( 1.0 - ( exp ( ( _v + 60.0 ) / 6.0 ) ) ) ;
		hBeta = ( - 0.015 * ( - _v - 60.0 ) ) / ( 1.0 - ( exp ( ( - _v - 60.0 ) / 6.0 ) ) );
		hInf = hAlpha / ( hAlpha + hBeta ) ;
        h=hInf;

    }
}

void output_states(NrnThread* _nt, int mech_id) {
    Mechanism _ml = _nt->ml[mech_id];
    int _cntml = _ml.nodecount;
    double * restrict _p = _ml.data;
    int _iml;

    printf("nodecount: %d\n", _cntml);

    for (_iml = 0; _iml < _cntml; ++_iml) {
        printf("NaTs2_t m: %g\n", m);
        printf("NaTs2_t h: %g\n", h);
    }
}

