function [ U1, Uend ] = computeBoundaries( wpIn, wpOut, wmIn, wmOut, param )
%UPDATEBOUNDARIES Summary of this function goes here
%   Detailed explanation goes here

    U1 = zeros(1,2);
    U1(1) = ( (wpIn - wmIn )./(8*sqrt( param.R0(1)) ) + (param.A0(1)).^(1/4) ).^4;
    U1(2) = ( (wpIn + wmIn) / 2 );
    
    Uend = zeros(1,2);
    Uend(1) = ( (wpOut - wmOut )./(8*sqrt( param.R0(end)) ) + (param.A0(end)).^(1/4) ).^4;
    Uend(2) = ( (wpOut + wmOut ) / 2 );
    

end

