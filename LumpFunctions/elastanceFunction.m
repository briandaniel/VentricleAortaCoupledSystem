function [ E ] = elastanceFunction( t, paramChamber )
%HILLFUNCTION Summary of this function goes here
%   Detailed explanation goes here

    E = paramChamber.k*hillFunctions( t, paramChamber ) + paramChamber.Emin;

end

