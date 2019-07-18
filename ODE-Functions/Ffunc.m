function [ F, Z ] = Ffunc( y, Z, wm, t, paramLA, paramMV, paramLV, param )
%FFUNC Summary of this function goes here
%   Detailed explanation goes here

    F = zeros(6,1);

    a1 = y(1);
    a2 = y(2);
    a3 = y(3);
    Pla = y(4);
    qmv = y(5);
    zetamv = y(6);
    
    % These values are required for the pressure/flow matching condition
    A0 = param.paramLVOT.A0(1);
    R0 = param.paramLVOT.R0(1);
    alphaLVOT = param.paramLVOT.alpha(1);
    
    [ alpha, kappa, eta, dVda, ~ ] = evaluateModel( a1, a2, a3, t, paramLV );
    
    % Newton iteration for the da/dt's + Plv
    k = 0;
    G = ones(4,1);
    while( k < 100 && norm(G,2) > 1e-10 )
        [ G ] = Gfunc( Z, alpha, kappa, eta, dVda, qmv, wm, A0, R0, alphaLVOT );
        [ J ] = Gjacobian( Z, alpha, eta, dVda, wm, A0, R0, alphaLVOT );

        Z = Z - J\G;
        k = k+1;
    end

    % Pull Z values
    F(1:3) = Z(1:3);
    Plv = Z(4);
    
    % dPla/dt
    F(4) = evaluateLA( Pla, qmv, param.Ppv, paramLA );

    % Valve time derivatives
  
    [ dqmv_dt, dzetamv_dt ] = ...
               computeValve( qmv, zetamv, Pla, Plv, paramMV, param );
	F(5) = dqmv_dt;
    F(6) = dzetamv_dt;

end


function [ G ] = Gfunc( Z, alpha, kappa, eta, dVda, qmv, wm, A0, R0, alphaLVOT )

    G = zeros(4,1);
    
    G(1) = alpha(1,1)*Z(1) + alpha(1,2)*Z(2) + alpha(1,3)*Z(3) - eta(1)*Z(4) + kappa(1);
    G(2) = alpha(2,1)*Z(1) + alpha(2,2)*Z(2) + alpha(2,3)*Z(3) - eta(2)*Z(4) + kappa(2);
    G(3) = alpha(3,1)*Z(1) + alpha(3,2)*Z(2) + alpha(3,3)*Z(3) - eta(3)*Z(4) + kappa(3);
    Ain = A0*( Z(4)/alphaLVOT + 1 )^2;
    G(4) = dVda(1)*Z(1) + dVda(2)*Z(2) - qmv ...
           + Ain*( wm + 4*sqrt(R0)*(Ain)^(1/4) - 4*sqrt(R0)*A0^(1/4) );

end


function [ J ] = Gjacobian( Z, alpha, eta, dVda, wm, A0, R0, alphaLVOT )

    J = zeros(4,4);
    
    J(1:3,1:3) = alpha;
    J(4,1:3) = dVda;
    J(1:3,4) = -eta;
    
    % This nonlinear derivative was checked with a quick finite difference check
    Ain = A0*( Z(4)/alphaLVOT + 1 )^2;   
    dG4dPlv = 2.*A0/alphaLVOT*( Z(4)/alphaLVOT + 1 )*( wm + 5*sqrt(R0)*(Ain)^(1/4) - 4*sqrt(R0)*A0^(1/4) );
    J(4,4) = dG4dPlv;
    

end












