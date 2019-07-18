function [ dV_da1, dV_da2 ] = volume_derivatives( subel,param)
%VOLUME_DERIVATIVES Summary of this function goes here
%   Detailed explanation goes here


% % Call backs
% a = subel.a;
% dmu_dnu  = subel.dmu_dnu ;
% dmu_da1  = subel.dmu_da1 ;
% dmu_da2 = subel.dmu_da2;
% ddmu_dnuda1 = subel.ddmu_dnuda1;
% ddmu_dnuda2  = subel.ddmu_dnuda2 ;
% mu = subel.mu;
% nu_vec = param.nu_vec; 
% 
% 
% % values only taken on inner wall
% mu = mu(1,:);
% nu = nu_vec;
% dmu_dnu = dmu_dnu(1,:);
% ddmu_dnuda1 = ddmu_dnuda1(1,:);
% ddmu_dnuda2 = ddmu_dnuda2(1,:);
% dmu_da2 = dmu_da2(1,:);
% dmu_da1 = dmu_da1(1,:);
% 
% s = sin(nu);
% s_sq = s.^2;
% sh = sinh(mu);
% c = cos(nu);
% sh_sq = sh.^2;
% 
% integrand1 = a^2/2*s_sq.*sh.*(-3.*s.*sinh(2*mu)+2*c.*sh_sq.*(3*dmu_dnu + a.*ddmu_dnuda1) ...
%                 - a*(s+3*cosh(2*mu).*s - 3.*c.*dmu_dnu.*sinh(2*mu) ).*dmu_da1);
% 
% integrand2 = a^3/2*s_sq.*sh.*(-2.*c.*sh_sq.*ddmu_dnuda2 + ...
%             (s + 3*cosh(2*mu).*s - 3*dmu_dnu.*c.*sinh(2*mu)).*dmu_da2);
%  
% dV_da1 = pi* numeric_integral(nu_vec,integrand1);
% dV_da2 = pi* numeric_integral(nu_vec,integrand2);



% Call backs
a = subel.a;
dmu_dnu  = subel.dmu_dnu ;
dmu_da1  = subel.dmu_da1 ;
dmu_da2 = subel.dmu_da2;
ddmu_dnuda1 = subel.ddmu_dnuda1;
ddmu_dnuda2  = subel.ddmu_dnuda2 ;
mu = subel.mu;
nu_vec = param.nu_vec; 

% values only taken on inner wall
mu_in = mu(1,:);



dmuin_da1 = subel.dmu_da1(1,:);

integrand = (1.0/12.0*cosh(3*mu_in) - 3.0/4.0*cosh(mu_in)).*sin(nu_vec) + ...
            cosh(mu_in).*sin(nu_vec).^3 - (-2.0/3.0*sin(nu_vec) + sin(nu_vec).^3);        
        
term1 = 6*pi*a.^2.*numeric_integral(nu_vec,integrand);

integrand_da1 = (1.0/12.0*sinh(3*mu_in)*3.*dmuin_da1 - 3.0/4.0*sinh(mu_in).*dmuin_da1).*sin(nu_vec) + ...
            sinh(mu_in).*dmuin_da1.*sin(nu_vec).^3;        
        
term2 = 2*pi*a.^3.*numeric_integral(nu_vec,integrand_da1);

dV_da1 = term1 + term2;


dmuin_da2 = subel.dmu_da2(1,:);

integrand_da2 = (1.0/4.0*sinh(3*mu_in) - 3.0/4.0*sinh(mu_in) + sinh(mu_in).*sin(nu_vec).^2  ).* dmuin_da2.* sin(nu_vec);
  
dV_da2 = 2*pi*a.^3.*numeric_integral(nu_vec,integrand_da2);


end

