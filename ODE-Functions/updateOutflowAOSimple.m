function [ Psa, wmOutAO ] = updateOutflowAOSimple( Psa, dt, U_AO, wpOutAO, param)
%UPDATEAOV Summary of this function goes here
%   Detailed explanation goes here
    
    % Pull the starting guess
    Aao = U_AO(end,1);
    uao = U_AO(end,2);
    
    Psa = Psa + dt/param.paramSA.C.*( Aao*uao - (Psa - param.Psp)/param.paramSA.R );

    X0 = Aao;
    
    % Evaluate the RHS at the previous time step
    Hnewton = @(X) Haov(X, Psa, wpOutAO, param );
    
    [ X ] = approximate_newton_iteration_Ndim( Hnewton, X0, param.newtParamAOV );

    
    % Compute incomming characteristic
    Aao = X;
    R0ao = param.paramAO.R0(end);
    A0ao = param.paramAO.A0(end);
    uao = wpOutAO - 4*R0ao^(1/2)*Aao^(1/4) + 4*R0ao^(1/2)*A0ao^(1/4);
    wmOutAO = uao - 4*R0ao^(1/2)*Aao^(1/4) + 4*R0ao^(1/2)*A0ao^(1/4);
    
    
end

function [ H ] = Haov( X, Psa_np1, wpOutAO, param )

        
    % downstream aortic pressure 
    A0 = param.paramAO.A0(end);
    R0 = param.paramAO.R0(end);
    Pao = param.paramAO.alpha(end)*( sqrt( X/A0) - 1 );
    Zc = param.paramAO.Zc_downstream;

    H = -Pao + Psa_np1 + Zc*(wpOutAO - 4*R0^(1/2)*X^(1/4) + 4*R0^(1/2)*A0^(1/4))*X;
    
end































          