function [ At] = activation_func( t, Ta, Tc, eps_fed, kp )
%ACTIVATION_FUNC Summary of this function goes here
%   Detailed explanation goes here

    % Reduced time
    tr = mod(t-(Tc-Ta),Tc);

    d = 1./(1+kp*eps_fed);

    At = 0*t;  
	At(tr <= Ta) = ( sin(pi*tr(tr <= Ta)/Ta) ).^d;
    

end


