function [ G, Vlv ] = computeDerivatives( a1, a2, a3, qmv, qlvot, t, paramLV )
%COMPUTEDERIVATIVES Summary of this function goes here
%   Detailed explanation goes here

    [ alpha, kappa, eta, dVda, Vlv ] = evaluateModel( a1, a2, a3, t, paramLV );
    
    A = zeros(4,4);
    A(1:3,1:3) = alpha;
    A(4,1:3) = dVda;
    A(1:3,4) = eta;
    
    b = zeros(4,1);
    b(1:3) = kappa;
    b(4) = qmv-qlvot;
    
    G = A\b;
    
end

