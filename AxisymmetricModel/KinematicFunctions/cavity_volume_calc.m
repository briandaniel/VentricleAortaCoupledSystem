function [ vol ] = cavity_volume_calc(  mu_in, nu_vec, dmu_dnu_in  , a)
%VOLUME_CALC Summary of this function goes here
%   Detailed explanation goes here

            
    integrand = (1.0/12.0*cosh(3*mu_in) - 3.0/4.0*cosh(mu_in)).*sin(nu_vec) + ...
            cosh(mu_in).*sin(nu_vec).^3 - (-2.0/3.0*sin(nu_vec) + sin(nu_vec).^3);

    vol = 2*pi*a.^3.*numeric_integral(nu_vec,integrand);

end

