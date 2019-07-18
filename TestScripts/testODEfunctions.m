
clear all
close all

run modelParameters.m
run LVparameters.m

%%
% Time domain
Nc = 1;
tmax = param.T*Nc;
dt = 0.0001;
t = 0;

% LVOT initial conditions
u0_LVOT = zeros(paramLVOT.Nx,1);
A0_LVOT = paramLVOT.A0;
U0_LVOT = [ A0_LVOT, u0_LVOT ];
U_LVOT = U0_LVOT;

% Initial values
a1 = 0;
a2 = 0;
a3 = 0;
Pla = 0;
qmv = 0; 
zetamv = 0;

% These should probably approximated before solution loop
Z = [0, 0, 0, 0]'; % Starting guess for derivatives, etc.

y0 = [ a1, a2, a3, Pla, qmv, zetamv ]';
y = y0;



%%
tic
[ y, U_LVOT, Z ] = ODEstep( t, y, U_LVOT, Z, dt, param );
toc



































































