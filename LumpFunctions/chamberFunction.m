function [ dVchamberdt, Pchamber ] = chamberFunction( treduced, Vchamber, qin, qout, paramChamber )
%CHAMBERFUNCTION Summary of this function goes here
%   Detailed explanation goes here

    qchamber = qin - qout;
    E = elastanceFunction( treduced, paramChamber);
    Pchamber = E*(Vchamber-paramChamber.V0)*(1+paramChamber.Ks*qchamber);
    dVchamberdt = qchamber;
    
end

