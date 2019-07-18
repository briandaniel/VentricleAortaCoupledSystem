function [ Se ] = elastic_stress ( C_prol, E_prol, sps, cps, paramLV)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
    [ C_fib ] = tensor_rotate_prolate_to_fiber( sps, cps, C_prol );
    [ E_fib ] = tensor_rotate_prolate_to_fiber( sps, cps, E_prol );


    bff = paramLV.bff;
    bxx = paramLV.bxx;
    bfx = paramLV.bfx;

    eW = exp( bff.*E_fib(:,:,3,3).^2 + bxx.*(E_fib(:,:,2,2).^2 ...
              + E_fib(:,:,1,1).^2 + 2.*E_fib(:,:,1,2).^2) ...
              + bfx.*(2*E_fib(:,:,3,1).^2 + 2.*E_fib(:,:,3,2).^2 ) );

    Se_fib = zeros(size(C_fib));     

    ke = paramLV.ke;
    Se_fib(:,:,1,1) = ke*eW.*bxx.*E_fib(:,:,1,1);
    Se_fib(:,:,2,2) = ke*eW.*bxx.*E_fib(:,:,2,2);
    Se_fib(:,:,3,3) = ke*eW.*bff.*E_fib(:,:,3,3);
    Se_fib(:,:,1,2) = ke*eW.*bxx.*E_fib(:,:,1,2);
    Se_fib(:,:,3,1) = ke*eW.*bfx.*E_fib(:,:,3,1);
    Se_fib(:,:,3,2) = ke*eW.*bfx.*E_fib(:,:,3,2);

    Se_fib(:,:,2,1) = Se_fib(:,:,1,2);
    Se_fib(:,:,1,3) = Se_fib(:,:,3,1);
    Se_fib(:,:,2,3) = Se_fib(:,:,3,2);

    % Rotate from fiber coordinates to prolate spheroidal coordinates
    Se = tensor_rotate_fiber_to_prolate(sps, cps, Se_fib);

end








