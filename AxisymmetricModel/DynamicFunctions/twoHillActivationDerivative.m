function [ AtPrime ] = twoHillActivationDerivative( t, m1, m2, tau1, tau2, Tc, Ts )
%TWOHILLACTIVATION Summary of this function goes here
%   Detailed explanation goes here

    t_cycle = mod(t - Ts, Tc);
    
    g1 = ( t_cycle / tau1).^m1;
    g2 = ( t_cycle / tau2).^m2;
    
    g1prime = (m1/tau1)*( t_cycle / tau1 ).^(m1-1); 
    g2prime = (m2/tau2)*( t_cycle / tau2 ).^(m2-1);
    
    AtPrime = g1prime./( (1+g1).^2.*(1+g2) ) - (g1*g2prime)./( (1+g1).*(1+g2).^2 );
    

end

