function [ U ] = vesselUpdateLW( U, dt, dx, paramVessel, x, xmid, t  )
%VESSELUPDATE Summary of this function goes here
%   Detailed explanation goes here


    F = flux(U, paramVessel.alpha, paramVessel.A0, paramVessel);
    S = source( U,paramVessel.mu,  paramVessel );
    
    Umid = 0.5*(U(2:end,:) + U(1:end-1,:)) - dt/(2*dx)*( F(2:end,:) - F(1:end-1,:) ) - ...
               dt/4*( S( 2:end,:) + S(1:end-1,:) );
    
    Fmid = flux( Umid, paramVessel.alpha_mid, paramVessel.A0_mid, paramVessel );
    smid = source( Umid, paramVessel.mu_mid, paramVessel );
    
    U(2:end-1,:) = U(2:end-1,:) - dt/dx.*( Fmid(2:end,:) - Fmid(1:end-1,:) ) - ...
        dt/2*( smid( 2:end,:) + smid(1:end-1,:) );

end



function [ F ] = flux( U, alpha, A0, paramVessel )

    F = zeros(size(U));
    F(:,1) = U(:,2).*U(:,1);
    F(:,2) = U(:,2).^2/2 + alpha/paramVessel.rho.*( sqrt(U(:,1)./A0) - 1 );

end


function [ C ] = source( U, mu, paramVessel ) 

    C = zeros(size(U));
    C(:,1) = 0;
    C(:,2) = paramVessel.xi*pi.*mu./paramVessel.rho .* ( U(:,2)./U(:,1) );

end

