function [ alpha, A0, R0, E, mu ] = alpha_A0_func( x, paramVessel )
%ALPHA_A0_FUNC Summary of this function goes here
%   Detailed explanation goes here


    A0 = pchip( paramVessel.xdata, paramVessel.A0data, x );
    h0 = pchip( paramVessel.xdata, paramVessel.h0data, x );
    E = pchip( paramVessel.xdata, paramVessel.Edata, x );
    mu = pchip( paramVessel.xdata, paramVessel.mudata, x );

%     plot(x,A0, paramVessel.xdata, paramVessel.A0data,'o' );

    % Put the spatial dependence of alpha/A0 here
    alpha = sqrt( pi ./ A0 ).*( E .*h0 )./( 1 - paramVessel.nu^2 );
    
    % Auxiliary
    R0 = alpha./(2*paramVessel.rho*A0.^(1/2));
    
    
end

