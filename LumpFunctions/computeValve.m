function [ dqvalve_dt, dZeta_dt ] = ...
                   computeValve( qvalve, zeta, upstreamP, downstreamP, paramValve, param )
%COMPUTEVALVE Summary of this function goes here
%   Detailed explanation goes here

    % Pressure difference is computed by subtracting the downstream
    % pressure from the upstream pressure
    deltaP = upstreamP - downstreamP;
    
    % Compute the effective areas
    Aeff = zeta*(paramValve.Aeff_max - paramValve.Aeff_min) + paramValve.Aeff_min;
    
    % Flow determining parameters
    B = param.rho/(2 * Aeff^2);
    L = param.rho.*paramValve.leff/(Aeff);
    
    % Flow derivative
    dqvalve_dt = 1/L*( deltaP - B*qvalve*abs(qvalve) );
    
    % Zeta derivative
    dZeta_dt = 0;
    if (deltaP > paramValve.deltaP_open)
        dZeta_dt = (1-zeta)*paramValve.Kvo*(deltaP -  paramValve.deltaP_open);
    elseif ( deltaP < paramValve.deltaP_close )
        dZeta_dt = zeta*paramValve.Kvc*(deltaP -  paramValve.deltaP_close);
    end
    
end

