function [ wpInLVOT ] = updateInflowLVOT( Plv, wmInLVOT, param )
%UPDATEINFLOWLVOT Summary of this function goes here
%   Detailed explanation goes here
    
    A0 = param.paramLVOT.A0(1);
    R0 = param.paramLVOT.R0(1);
    Alvot = param.paramLVOT.A0(1)*( Plv/param.paramLVOT.alpha(1) + 1 )^2;
    ulvot = wmInLVOT + 4*sqrt(R0)*(Alvot)^(1/4) - 4*sqrt(R0)*A0^(1/4);
    wpInLVOT = ulvot + 4*R0.^(1/2).*Alvot^(1/4) - 4*R0.^(1/2)*A0^(1/4);
    

end

