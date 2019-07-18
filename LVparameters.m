run modelParameters.m

addpath('AxisymmetricModel')
addpath('AxisymmetricModel\KinematicFunctions')
addpath('AxisymmetricModel\DynamicFunctions')
addpath('AxisymmetricModel/KinematicFunctions')
addpath('AxisymmetricModel/DynamicFunctions')



%%
% Physical domain parameters
paramLV.Nmu = 10;                               % Number of points in the \mu direction
paramLV.Nnu = 10;                               % Number of points in the \nu direction
paramLV.nu_up = 1.4;                            % Smallest \nu value ( \nu ranges in [ nu_up, pi] )
paramLV.nu_vec = linspace(paramLV.nu_up,pi,paramLV.Nnu);        % \nu vector

% Ventricular contraction
paramLV.Tc = param.T;                         	% Length of full cardiac cycle (s) 
paramLV.Ta = 0.5;                        	% Length of active contraction (s)

% Initial structure parameters
paramLV.a0 = 5.08;                              % Initial spheroid length (cm)
paramLV.muin0 = .37; %matched                   % Inner wall starting \mu value
paramLV.muout0 = .6; %matched                 % Outer wall starting \mu value

% Muscle fiber parametrization angles
paramLV.psi_in_b0 = 85 * (pi / 180); %rad       % Angle of fiber wrapping on endocardial surface 
paramLV.psi_out_b0 = -65 * (pi / 180); %rad      % Angle of fiber wrapping on epidcardial surface

% Active stress parameters
% paramLV.km = 60e4; % dynes/cm^2                  % Initial strength of active fiber stress
paramLV.km = 60e4; % dynes/cm^2                  % Initial strength of active fiber stress

paramLV.kav = 5e4; % dyne/cm^2 * s               % Coefficient of the viscous component of the active fiber stress
paramLV.kp = 0;                                 % Coefficient of the end-diastolic response (pre-load dependent active stress)
paramLV.ked = 0.0;                              % Coefficient of the end-diastolic response (pre-load dependent active stress)
paramLV.m = 0.3;                                % Slope of the stress/strain relationship in the active fibers

% Viscosity
paramLV.kv = 0.1 * 10^4; % Converted from kPa*s to dyne/cm^2 * s

% Elasticity
paramLV.ke = 1.5 * 10^4; % Converted from kPa to dyne/cm^2
paramLV.bff = 2.5;
paramLV.bfx = 2.5;
paramLV.bxx = 2.5;

% Assign to main parameter struct
param.paramLV = paramLV;








