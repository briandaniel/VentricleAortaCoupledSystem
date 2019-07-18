function [ wm ] = updateOutflowLVOT( Alvot, wp, qaov, param )
%UPDATEOUTFLOWLVOT Summary of this function goes here
%   Detailed explanation goes here

    %  Use a newton iteration to solve for Alvot.
    R0 = param.paramLVOT.R0(end);
    A0 = param.paramLVOT.A0(end);
    k = 1;
    F = 1;
    while ( k < 100 && norm(F) > 1e-10 )
        
       [ F, dFdA ] = newtonFunction( Alvot, wp, qaov, R0, A0 );
       Alvot = Alvot - F/dFdA;

       k = k+1; 

    end
    
    ulvot = wp - 4*R0^(1/2)*Alvot^(1/4) + 4*R0^(1/2)*A0^(1/4);
    wm = ulvot - 4*R0^(1/2)*Alvot^(1/4) + 4*R0^(1/2)*A0^(1/4);
    
end




function [ F, dFdA ] = newtonFunction( Alvot, wp, qaov, R0, A0 )

    ulvot = wp - 4*R0^(1/2)*Alvot^(1/4) + 4*R0^(1/2)*A0^(1/4);
    F = -qaov + Alvot*(ulvot);
    dFdA = ulvot - Alvot^(1/4)*R0^(1/2);
    
end