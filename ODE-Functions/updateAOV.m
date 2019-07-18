function [ w, wmOutLVOT, wpInAO ] = updateAOV( wn, dt, U_LVOT, U_AO, wpOutLVOT, wmInAO, param)
%UPDATEAOV Summary of this function goes here
%   Detailed explanation goes here
    
    % Pull the starting guess
    qaov = wn(1);
    zetaaov = wn(2);
    Alvot = U_LVOT(end,1);
    Aao = U_AO(1,1);
    
    X0 = [ qaov, zetaaov, Alvot, Aao ]';
    
    % Evaluate the RHS at the previous time step
    [ Fn ] = Faov( X0, wpOutLVOT, wmInAO, param );

    Hnewton = @(X) Haov(X, wn, Fn, dt, wpOutLVOT, wmInAO, param );
    
    [ X ] = approximate_newton_iteration_Ndim( Hnewton, X0, param.newtParamAOV );
    
    w = wn;
    w(1:2) = X(1:2);
    Alvot = X(3);
    Aao = X(4);
    
    % Compute characteristics
    R0lvot = param.paramLVOT.R0(end);
    A0lvot = param.paramLVOT.A0(end);
    ulvot = wpOutLVOT - 4*R0lvot^(1/2)*Alvot^(1/4) + 4*R0lvot^(1/2)*A0lvot^(1/4);
    wmOutLVOT = ulvot - 4*R0lvot^(1/2)*Alvot^(1/4) + 4*R0lvot^(1/2)*A0lvot^(1/4);
    
    R0ao = param.paramAO.R0(1);
    A0ao = param.paramAO.A0(1);
    uao = wmInAO + 4*R0ao^(1/2)*Aao^(1/4) - 4*R0ao^(1/2)*A0ao^(1/4);
    wpInAO = uao + 4*R0ao^(1/2)*Aao^(1/4) - 4*R0ao^(1/2)*A0ao^(1/4);

    
end

function [ H ] = Haov( X, wn, Fn, dt, wpLVOT, wmAO, param )

    [ F ] = Faov( X, wpLVOT, wmAO, param );
   
    H = zeros(4,1);
    H(1:2) = - X(1:2) + wn(1:2) + dt/2*( F(1:2) + Fn(1:2) );
    H(3:4) = F(3:4);

end


function [ F ] = Faov( X, wpLVOT, wmAO, param )

    qaov = X(1);
    zetaaov = X(2);
    Alvot = X(3);
    Aao = X(4);
    
    Plvot = param.paramLVOT.alpha(end)*( sqrt( Alvot/param.paramLVOT.A0(end)) - 1 );
    Pao = param.paramAO.alpha(1)*( sqrt( Aao/param.paramAO.A0(1)) - 1 );
    
    F = zeros(4,1);
    [ dqaov_dt, dzetaaov_dt ] = ...
               computeValve( qaov, zetaaov, Plvot, Pao, param.paramAOV, param );
	F(1) = dqaov_dt;
    F(2) = dzetaaov_dt;
   
    R0lvot = param.paramLVOT.R0(end);
    A0lvot = param.paramLVOT.A0(end);
    F(3) =  -qaov + Alvot*( wpLVOT - 4*sqrt(R0lvot)*Alvot^(1/4) + 4*sqrt(R0lvot)*A0lvot^(1/4) );
    
      
    R0ao = param.paramAO.R0(1);
    A0ao = param.paramAO.A0(1);
    F(4) =  -qaov + Aao*( wmAO + 4*sqrt(R0ao)*Aao^(1/4) - 4*sqrt(R0ao)*A0ao^(1/4) );
    
end
































          