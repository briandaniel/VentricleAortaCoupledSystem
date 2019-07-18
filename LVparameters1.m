
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
paramLV.eps_fed = 0.0; % This can/should be adjusted

% Viscosity
paramLV.kv = 0.1 * 10^4; % Converted from kPa*s to dyne/cm^2 * s

% Elasticity
paramLV.ke = 1.5 * 10^4; % Converted from kPa to dyne/cm^2
paramLV.bff = 2.5;
paramLV.bfx = 2.5;
paramLV.bxx = 2.5;


% Alternative activation (two-hill curve)
paramLV.useTwoHill = 1; % Set to 1 if you want to use the activation function
paramLV.m1 = 5.32;
paramLV.m2 = 8;
paramLV.tau1 = 0.10;
paramLV.tau2 = 0.3;
paramLV.Ts = 0.4;
[tMax,paramLV.hillMaxVal] = computeHillMax( paramLV.m1, paramLV.m2, paramLV.tau1, paramLV.tau2, paramLV.Tc);

% Plot the activation curves
t = linspace(0,2,1e3);
[ At] = activation_func( t, paramLV.Ta, paramLV.Tc, 0, paramLV.kp );
At2 = twoHillActivation( t, paramLV.m1, paramLV.m2, paramLV.tau1, paramLV.tau2, paramLV.Tc, paramLV.Ts, paramLV.hillMaxVal );
plot(t,At,t,At2);


% Assign to main parameter struct
param.paramLV = paramLV;



