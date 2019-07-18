function [ Psa, wmOutAO ] = updateOutflowAO( Psa, dt, U_AO, wpOutAO, param)
%UPDATEAOV Summary of this function goes here
%   Detailed explanation goes here
    
    % Pull the starting guess
    Aao = U_AO(end,1);
    
    X0 = [ Psa, Aao ]';
    
    % Evaluate the RHS at the previous time step
    [ Fn ] = Fsa( X0, param );

    Hnewton = @(X) Haov(X, Psa, Fn, dt, wpOutAO, param );
    
    [ X ] = approximate_newton_iteration_Ndim( Hnewton, X0, param.newtParamAOV );
    
    % New value of Psa
    Psa = X(1);
    
    % Compute incomming characteristic
    Aao = X(2);
    R0ao = param.paramAO.R0(end);
    A0ao = param.paramAO.A0(end);
    uao = wpOutAO - 4*R0ao^(1/2)*Aao^(1/4) + 4*R0ao^(1/2)*A0ao^(1/4);
    wmOutAO = uao - 4*R0ao^(1/2)*Aao^(1/4) + 4*R0ao^(1/2)*A0ao^(1/4);
    
    
end

function [ H ] = Haov( X, Psa_n, Fn, dt, wpOutAO, param )

    [ F ] = Fsa( X, param );
   
    H = zeros(2,1);
    H(1) = - X(1) + Psa_n + dt/2*( F + Fn );
    
    
        
    % downstream aortic pressure 
    A0 = param.paramAO.A0(end);
    R0 = param.paramAO.R0(end);
    Pao = param.paramAO.alpha(end)*( sqrt( X(2)/A0) - 1 );
    Zc = param.paramAO.Zc_downstream;

    H(2) = -Pao + X(1) + Zc*(wpOutAO - 4*R0^(1/2)*X(2)^(1/4) + 4*R0^(1/2)*A0^(1/4))*X(2);
    
end


function [ F ] = Fsa( X, param )

    Psa = X(1);
    Aao = X(2);
    
    % downstream aortic pressure 
    A0 = param.paramAO.A0(end);
    Pao = param.paramAO.alpha(end)*( sqrt( Aao/A0) - 1 );

    % RHS of lump equation for Psa
    Zc = param.paramAO.Zc_downstream;
    F = 1/param.paramSA.C*( ( Pao - Psa )/Zc - (Psa - param.Psp)/param.paramSA.R );

end
































          