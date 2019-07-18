
close all
clear all

run LVparameters.m;

%%

t = 0;
a1 = 0.1;
a2 = 0.2;
a3 = -0.1;

[ alpha, kappa, eta, dVda, Vlv ] = evaluateModel( a1, a2, a3, t, paramLV )


qmv = 0.1;
qlvot = 0.1;
[ G ] = computeDerivatives( a1, a2, a3, qmv, qlvot, t, paramLV )





























