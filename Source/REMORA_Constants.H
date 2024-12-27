#ifndef REMORA_CONSTANTS_H_
#define REMORA_CONSTANTS_H_

#include <AMReX_REAL.H>

constexpr amrex::Real PI = amrex::Real(3.14159265358979323846264338327950288);
//Note that amrex source code uses fewer digits of pi:
//     constexpr Real pi = Real(3.1415926535897932);
//TODO: Change these types of macros to 'const'
constexpr amrex::Real R_d      = amrex::Real(287.0);   // gas constant for dry air [J/(kg-K)]
constexpr amrex::Real R_v      = amrex::Real(461.6);   // gas constant for dry air [J/(kg-K)]
constexpr amrex::Real c_p      = amrex::Real(1004.5);   // specific heat at constant pressure for dry air [J/(kg-K)]
constexpr amrex::Real p_0        = amrex::Real(1.0e5); // reference surface pressure [Pa]
constexpr amrex::Real Gamma      = amrex::Real(1.4);   // C_p/C_v [-]
constexpr amrex::Real KAPPA      = amrex::Real(0.41);  // von Karman constant

// Derived Constants
constexpr amrex::Real ip_0    = amrex::Real(1./p_0);
constexpr amrex::Real iR_d    = amrex::Real(1./R_d);
constexpr amrex::Real iGamma  = amrex::Real(1./Gamma);
constexpr amrex::Real rdOcp   = amrex::Real(R_d/c_p);

//  Polynomial  expansion  coefficients for the computation of
//  "in situ" density  and other associated quantities via the
//  nonlinear equation of state for seawater  as a function of
//  potential temperature, salinity, and pressure (Jackett and
//  McDougall, 1992).

constexpr amrex::Real A00 = amrex::Real(1.909256e+04);
constexpr amrex::Real A01 = amrex::Real(2.098925e+02);
constexpr amrex::Real A02 = amrex::Real(-3.041638e+00);
constexpr amrex::Real A03 = amrex::Real(-1.852732e-03);
constexpr amrex::Real A04 = amrex::Real(-1.361629e-05);
constexpr amrex::Real B00 = amrex::Real(1.044077e+02);
constexpr amrex::Real B01 = amrex::Real(-6.500517e+00);
constexpr amrex::Real B02 = amrex::Real(1.553190e-01);
constexpr amrex::Real B03 = amrex::Real(2.326469e-04);
constexpr amrex::Real D00 = amrex::Real(-5.587545e+00);
constexpr amrex::Real D01 = amrex::Real(7.390729e-01);
constexpr amrex::Real D02 = amrex::Real(-1.909078e-02);
constexpr amrex::Real E00 = amrex::Real(4.721788e-01);
constexpr amrex::Real E01 = amrex::Real(1.028859e-02);
constexpr amrex::Real E02 = amrex::Real(-2.512549e-04);
constexpr amrex::Real E03 = amrex::Real(-5.939910e-07);
constexpr amrex::Real F00 = amrex::Real(-1.571896e-02);
constexpr amrex::Real F01 = amrex::Real(-2.598241e-04);
constexpr amrex::Real F02 = amrex::Real(7.267926e-06);
constexpr amrex::Real G00 = amrex::Real(2.042967e-03);
constexpr amrex::Real G01 = amrex::Real(1.045941e-05);
constexpr amrex::Real G02 = amrex::Real(-5.782165e-10);
constexpr amrex::Real G03 = amrex::Real(1.296821e-07);
constexpr amrex::Real H00 = amrex::Real(-2.595994e-07);
constexpr amrex::Real H01 = amrex::Real(-1.248266e-09);
constexpr amrex::Real H02 = amrex::Real(-3.508914e-09);
constexpr amrex::Real Q00 = amrex::Real(9.99842594e+02);
constexpr amrex::Real Q01 = amrex::Real(6.793952e-02);
constexpr amrex::Real Q02 = amrex::Real(-9.095290e-03);
constexpr amrex::Real Q03 = amrex::Real(1.001685e-04);
constexpr amrex::Real Q04 = amrex::Real(-1.120083e-06);
constexpr amrex::Real Q05 = amrex::Real(6.536332e-09);
constexpr amrex::Real U00 = amrex::Real(8.24493e-01);
constexpr amrex::Real U01 = amrex::Real(-4.08990e-03);
constexpr amrex::Real U02 = amrex::Real(7.64380e-05);
constexpr amrex::Real U03 = amrex::Real(-8.24670e-07);
constexpr amrex::Real U04 = amrex::Real(5.38750e-09);
constexpr amrex::Real V00 = amrex::Real(-5.72466e-03);
constexpr amrex::Real V01 = amrex::Real(1.02270e-04);
constexpr amrex::Real V02 = amrex::Real(-1.65460e-06);
constexpr amrex::Real W00 = amrex::Real(4.8314e-04);

#endif
