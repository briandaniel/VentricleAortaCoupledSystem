clear all
close all

addpath('LumpFunctions')
addpath('UtilityFunctions')
addpath('ODE-Functions')
addpath('Vessel-1D-functions')
addpath('TestScripts')


%% 
% Note: 1 dyne/cm^2 = 1 g/(cm*s^2), which are the actual units of pressure
% computed
% Unit conversion is 7.5e-4 mmHg = 1 dyne/cm^2 = 1 g/(cm*s^2)
% or 1 mmHg = 133 Pa = 1333.3 g/(cm*s^2) = 1333.3 dyne/cm^2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nonspecific %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.T = 0.89; % cycle duration s
param.rho = 1.06; % gm/cm^3 = gm/mL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Aortic Valve parameters %%%%%%%%%%%%%%%%%%%%%%%%%
paramAOV.Aeff_min = 0.0001;
paramAOV.Aeff_max = 6; % cm^2
paramAOV.leff = 0.2; %cm
% paramAOV.Kvo = 0.012; % 1/( s * dyne/cm^2 ) 
% paramAOV.Kvc = 0.012; %  1/( s * dyne/cm^2 ) 

paramAOV.Kvo = 0.15; % 1/( s * dyne/cm^2 ) 
paramAOV.Kvc = 0.15; %  1/( s * dyne/cm^2 ) 
paramAOV.deltaP_open = 0;
paramAOV.deltaP_close = 0;
param.paramAOV = paramAOV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% Mitral Valve parameters %%%%%%%%%%%%%%%%%%%%%%%%%
paramMV.Aeff_min = 0.0001; % cm^2
paramMV.Aeff_max = 6; % cm^2
paramMV.leff = 0.2; %cm
paramMV.Kvo = 0.15 ; % 1/( s * dyne/cm^2 ) 
paramMV.Kvc = 0.15 ; % 1/( s * dyne/cm^2 ) 
paramMV.deltaP_open = 0;
paramMV.deltaP_close = 0;

param.paramMV = paramMV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LVOT parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramLVOT.xbnd = [0,2]; % cm
paramLVOT.dx = 1; % cm
paramLVOT.x = ( paramLVOT.xbnd(1):paramLVOT.dx:paramLVOT.xbnd(2) )';
paramLVOT.xmid = ( ( paramLVOT.xbnd(1) + paramLVOT.dx/2 ):paramLVOT.dx:(paramLVOT.xbnd(2) - paramLVOT.dx/2) )';
paramLVOT.Nx = length(paramLVOT.x);
paramLVOT.rho = 1.06; % g/cm^3 = g/mL
paramLVOT.nu = 1/2; % For an incompressible fluid nu = 1/2
paramLVOT.xi = 0;
paramLVOT.inletmu = 0;
paramLVOT.E = 4*10^6; % g/(cm*s^2) = dyne/cm^2
paramLVOT.inletA0 = 4.5; % cm^2
paramLVOT.inleth0 = 0.16; % cm

[ paramLVOT.alpha, paramLVOT.A0, paramLVOT.R0, paramLVOT.mu  ] = alpha_A0_func_constant( paramLVOT.x, paramLVOT );
[ paramLVOT.alpha_mid, paramLVOT.A0_mid, paramLVOT.R0_mid, paramLVOT.mu_mid  ] = alpha_A0_func_constant( paramLVOT.xmid, paramLVOT );
param.paramLVOT = paramLVOT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AO parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramAO.xbnd = [0,45]; % cm
paramAO.dx = 1; % cm
paramAO.x = ( paramAO.xbnd(1):paramAO.dx:paramAO.xbnd(2) )';
paramAO.xmid = ( ( paramAO.xbnd(1) + paramAO.dx/2 ):paramAO.dx:(paramAO.xbnd(2) - paramAO.dx/2) )';
paramAO.Nx = length(paramAO.x);
paramAO.rho = 1.06; % g/cm^3 = g/mL
paramAO.nu = 1/2; % For an incompressible fluid nu = 1/2
paramAO.xi = 8;




paramAO.xdata = [5, 40];
paramAO.Edata = [1, 1]*2*10^6;
paramAO.A0data = [6.5, 3];
paramAO.h0data = [1, 1]*0.16;
paramAO.mudata = [1, 1]*0.04;



paramAO.coarct_pos = 20;
paramAO.coarct_width = 10;
paramAO.coarct_area_change =0;
paramAO.coarct_E_change = 0*10^6;


[ paramAO.alpha, paramAO.A0, paramAO.R0, paramAO.E, paramAO.mu ] = alpha_A0_func_linear_coarct( paramAO.x, paramAO );
% [ paramAO.alpha, paramAO.A0, paramAO.R0, paramAO.E, paramAO.mu ] = alpha_A0_func_linear( paramAO.x, paramAO );
[ paramAO.alpha_mid, paramAO.A0_mid, paramAO.R0_mid, paramAO.E_mid, paramAO.mu_mid ] = alpha_A0_func_linear_coarct( paramAO.xmid, paramAO );

figure('outerposition',[1500,100,600,1100]);
subplot(4,1,1)
hold on
plot(paramAO.xmid, paramAO.A0_mid)
plot(paramAO.x, paramAO.A0)
title('Area')

subplot(4,1,2)
hold on
plot(paramAO.xmid, paramAO.E_mid)
plot(paramAO.x, paramAO.E)
title('Elasticity')

subplot(4,1,3)
hold on
plot(paramAO.xmid, paramAO.alpha_mid)
plot(paramAO.x, paramAO.alpha)
title('alpha')

subplot(4,1,4)
hold on
plot(paramAO.xmid, paramAO.R0_mid)
plot(paramAO.x, paramAO.R0)
title('R0')


paramAO.Zc_downstream = sqrt( paramAO.alpha(end) * param.rho )/(paramAO.A0(end)*sqrt(2));
param.paramAO = paramAO;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lump models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upstream/downstream conditions
param.Ppv = 10 * 1333.3; % mmHg converted to dyne/cm^2
param.Psp = 50 * 1333.3; % mmHg converted to dyne/cm^2

% Left atrium windkessel
paramLA.C = 4.5e-3; % (mL) / (0.1Pa)
paramLA.R = 5; % (0.1Pa) / (mL/s) 
param.paramLA = paramLA;

% Systemic artery windkessel
paramSA.C = 0.01;
paramSA.R = 550; % (dyne/cm^2) / (mL/s) 

% paramSA.C = 1e-3;
% paramSA.R = 5; % (dyne/cm^2) / (mL/s) 

param.paramSA = paramSA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% Newton solver parameters %%%%%%%%%%%%%%%%%%%%%%%%%
newtParam.dx = 1e-7;
newtParam.xMinDiff = 1e-5;
newtParam.fNewtMin = 1e-5;
newtParam.maxIter = 1;
param.newtParam= newtParam;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% Newton solver parameters %%%%%%%%%%%%%%%%%%%%%%%%%
newtParamAOV.dx = 1e-7;
newtParamAOV.xMinDiff = 1e-5;
newtParamAOV.fNewtMin = 1e-5;
newtParamAOV.maxIter = 10;
param.newtParamAOV = newtParamAOV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


















