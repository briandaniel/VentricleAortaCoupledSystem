function [ alpha, A0, R0, mu ] = alpha_A0_func_constant( x, paramVessel )
%ALPHA_A0_FUNC Summary of this function goes here
%   Detailed explanation goes here

    A0 = paramVessel.inletA0*ones(size(x));
    h0 = paramVessel.inleth0*ones(size(x));
    mu = paramVessel.inletmu*ones(size(x));

    % Put the spatial dependence of alpha/A0 here
    alpha = sqrt( pi ./ A0 ).*( paramVessel.E .*h0 ./( 1 - paramVessel.nu^2 ) );
    
    % Auxiliary
    R0 = alpha./(2*paramVessel.rho*A0.^(1/2));


end

