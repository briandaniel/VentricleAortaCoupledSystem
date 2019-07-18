function [ H, Z ] = Hfunc( y, Z, t, dt, Fn, yn, param, wmIn )
%HFUNC Summary of this function goes here
%   Detailed explanation goes here

    [ Fnp1, Z ] = Ffunc( y, Z, wmIn, t, param.paramLA, param.paramMV, param.paramLV, param );

    % H is simply the finite difference approximation
    H = -y + dt/2*Fnp1 + yn + dt/2*Fn;
    
end

