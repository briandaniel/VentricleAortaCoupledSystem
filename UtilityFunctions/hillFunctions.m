function [ H ] = hillFunctions( t, paramChamber )
%HILLFUNCTIONS Summary of this function goes here
%   Detailed explanation goes here

    g1 = ( t/paramChamber.tau1 ).^paramChamber.m1;
    g2 = ( t/paramChamber.tau2 ).^paramChamber.m2;
    
    H = (g1./(1+g1)).*(1./(1 + g2));


end

