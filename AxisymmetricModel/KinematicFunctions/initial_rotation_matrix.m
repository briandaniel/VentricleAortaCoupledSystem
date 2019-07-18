function [ cps, sps, omega ] ... 
    = initial_rotation_matrix( subel,param )
%ROTATION_MATRIX Summary of this function goes here
%   Detailed explanation goes here

% Call backs
sh0 = subel.sh0;
ch0 = subel.ch0;
s = subel.s;
c = subel.c;
Nmu = param.Nmu;
Nnu = param.Nnu;
psi_in_b0 = param.psi_in_b0;
psi_out_b0 = param.psi_out_b0;
muin0 = param.muin0;
muout0 = param.muout0;
muvec0 = subel.muvec0;

% Calculate the initial psi angle values as a function of initial mu values
psi_eq = 1./(muout0 - muin0).*( psi_in_b0.*(muout0 - muvec0) ...
                                  + psi_out_b0.*(muvec0 - muin0) );
psiStar_eq = psi_eq;
psiStar_eq(psiStar_eq < 0) = pi + psi_eq(psiStar_eq < 0);

omega = 1./(tan(psiStar_eq).*tanh(muvec0));
omg = repmat( omega', 1, Nnu );

% Calculate psi
cps = sh0.*s.^2.*omg./ ...
    sqrt( sh0.^2 + s.^2 + sh0.^2.*s.^4.*omg.^2 );   

sps = sqrt( sh0.^2 + s.^2 )./ ...
    sqrt( sh0.^2 + s.^2 + sh0.^2.*s.^4.*omg.^2 );       

end


