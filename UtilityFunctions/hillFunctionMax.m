function [ k ] = hillFunctionMax( paramChamber, param )
%HILLFUNCTIONMAX Summary of this function goes here
%   Detailed explanation goes here


    obj = @(t) -hillFunctions( t, paramChamber );
    
    tmax = fminbnd(obj,0,param.T);
    
    
    k = (paramChamber.Emax - paramChamber.Emin)./hillFunctions( tmax, paramChamber );
    
    
    

end

