function [ Sf_const_prol ] = active_fiber_stress_constant( At, km, ked, m, ...
                        ell, ell0, eps_ed, sps, cps, F )
%ACTIVE_FIBER_STRESS Summary of this function goes here
%   Detailed explanation goes here


eps_f = 1/2*(ell.^2./ell0.^2 - 1 );

% active fiber stress component in the ff direction (All other directions
% are zero)
sigmaF = At.*km.*(1+m.*eps_f).*(1+ked.*eps_ed);

% Convert to Piola-Kirchoff stress tensor
SF = ell0.^2./ell.^2 .*sigmaF;

% Rotate from fiber coordinates back to prolate spheroidal coordaintes
[Nmu,Nnu] = size(ell);
Sf_const_fib = zeros(Nmu,Nnu,3,3);
Sf_const_fib(:,:,3,3) = SF;

Sf_const_prol = tensor_rotate_fiber_to_prolate(sps, cps, Sf_const_fib);

end



























