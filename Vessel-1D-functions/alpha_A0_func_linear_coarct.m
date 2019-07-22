function [ alpha, A0, R0, E, mu ] = alpha_A0_func_linear_coarct( x, paramVessel )
%ALPHA_A0_FUNC Summary of this function goes here
%   Detailed explanation goes here

    [A0, E] = computeCoarctValues(x,paramVessel);
    
    h0 = computeVariableValues( paramVessel.xdata, paramVessel.h0data, x );
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


function [A0, E] = computeCoarctValues( x, paramVessel )


    A0 = computeVariableValues( paramVessel.xdata, paramVessel.A0data, x );
    E = computeVariableValues( paramVessel.xdata, paramVessel.Edata, x );
    
    x_coarct = [x(1), paramVessel.coarct_pos + [ -paramVessel.coarct_width, - paramVessel.coarct_width*.01,...
                0,  paramVessel.coarct_width*.01, paramVessel.coarct_width ], x(end)];
    delA0_x_coarct = [0, 0,  paramVessel.coarct_area_change, ...
        paramVessel.coarct_area_change, paramVessel.coarct_area_change, 0, 0];
    delE_x_coarct = [0, 0,  paramVessel.coarct_E_change, ...
        paramVessel.coarct_E_change, paramVessel.coarct_E_change, 0, 0];
    
    
    delA0 = pchip(x_coarct,delA0_x_coarct,x);
    delE = pchip(x_coarct,delE_x_coarct,x);

    A0 = A0 - delA0;
    E = E + delE;
    
      
end













