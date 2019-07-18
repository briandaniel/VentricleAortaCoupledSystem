function [ vol ] = cavity_volume_calc_minimum_input( a1, a2, paramLV )
%VOLUME_CALC Summary of this function goes here
%   Detailed explanation goes here

    nu_vec = linspace(paramLV.nu_up,pi,20);
    mu_in = muCalc( paramLV.a0, a1, a2, paramLV.muin0, nu_vec, paramLV.muin0);
    
    a = paramLV.a0+a1;
            
    integrand = (1.0/12.0*cosh(3*mu_in) - 3.0/4.0*cosh(mu_in)).*sin(nu_vec) + ...
            cosh(mu_in).*sin(nu_vec).^3 - (-2.0/3.0*sin(nu_vec) + sin(nu_vec).^3);

    vol = 2*pi*a.^3.*numeric_integral(nu_vec,integrand);

end

