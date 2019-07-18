function [ wp ] = updateInflowAO( Aao, wm, qaov, param )
%UPDATEOUTFLOWLVOT Summary of this function goes here
%   Detailed explanation goes here

    %  Use a newton iteration to solve for Alvot.
    R0 = param.paramAO.R0(1);
    A0 = param.paramAO.A0(1);
    k = 1;
    F = 1;
    while ( k < 100 && norm(F) > 1e-10 )
        
       [ F, dFdA ] = newtonFunction( Aao, wm, qaov, R0, A0 );
       Aao = Aao - F/dFdA;

       k = k+1; 

    end
    
    uao = wm + 4*R0^(1/2)*Aao^(1/4) - 4*R0^(1/2)*A0^(1/4);
    wp = uao + 4*R0^(1/2)*Aao^(1/4) - 4*R0^(1/2)*A0^(1/4);
    
end




function [ F, dFdA ] = newtonFunction( Aao, wm, qaov, R0, A0 )

    uao = wm + 4*R0^(1/2)*Aao^(1/4) - 4*R0^(1/2)*A0^(1/4);
    F = -qaov + Aao*(uao);
    dFdA = uao + Aao^(1/4)*R0^(1/2);
    
end