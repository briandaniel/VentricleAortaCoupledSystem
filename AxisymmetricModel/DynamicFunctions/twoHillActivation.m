function [ At ] = twoHillActivation( t, m1, m2, tau1, tau2, Tc, Ts, hillMaxVal )
%TWOHILLACTIVATION Summary of this function goes here
%   Detailed explanation goes here

    t_cycle = mod(t - Ts, Tc);
    
    g1 = ( t_cycle / tau1).^m1;
    g2 = ( t_cycle / tau2).^m2;
    At = ( g1./(1 + g1).*(1./(1+g2)) )/hillMaxVal;
    

end

