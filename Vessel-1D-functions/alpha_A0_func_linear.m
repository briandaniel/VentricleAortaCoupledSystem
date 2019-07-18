function [ alpha, A0, R0, E, mu ] = alpha_A0_func_linear( x, paramVessel )
%ALPHA_A0_FUNC Summary of this function goes here
%   Detailed explanation goes here

    A0 = computeVariableValues( paramVessel.xdata, paramVessel.A0data, x );
    h0 = computeVariableValues( paramVessel.xdata, paramVessel.h0data, x );
    E = computeVariableValues( paramVessel.xdata, paramVessel.Edata, x );
    mu = computeVariableValues( paramVessel.xdata, paramVessel.mudata, x );


    % Put the spatial dependence of alpha/A0 here
    alpha = sqrt( pi ./ A0 ).*( E .*h0 )./( 1 - paramVessel.nu^2 );
    
    % Auxiliary
    R0 = alpha./(2*paramVessel.rho*A0.^(1/2));
    
    
end

function [f] = computeVariableValues( xData, fData, x )

    
    % regular domain 
    ux = ( x - xData(1) )./(xData(2) - xData(1) );
    f = (1 - ux).*fData(1) + ux.*fData(2);
    
    % adjust constant bounadries
    f( x <= xData(1) ) = fData(1);
    f( x >= xData(2) ) = fData(2);
%     keyboard

end
