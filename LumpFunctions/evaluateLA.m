function [ dPladt ] = evaluateLA( Pla, qmv, Ppv, paramLA )
%EVALUATELA Summary of this function goes here
%   Detailed explanation goes here

    dPladt = 1/paramLA.C * ( (Ppv - Pla)/paramLA.R - qmv );

end

