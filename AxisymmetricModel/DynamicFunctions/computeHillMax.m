function [tMax,hillMax] = computeHillMax( m1, m2, tau1, tau2, Tc)
%COMPUTEHILLMAX Summary of this function goes here
%   Detailed explanation goes here
 
    Ts = 0;
    
    % hill function has inflections, so safer to simply find good starting
    % bounds
    Nguess = 50;
    tMax = 0;
    hillMax = -1e10;
    for k=1:Nguess
        t = Tc*(k/Nguess);
        hillk = twoHillActivation( t, m1, m2, tau1, tau2, Tc, Ts, 1 );
        if(hillk > hillMax)
            hillMax = hillk;
            tMax = k/Nguess*Tc;
        end
    end
    
    
    ta = tMax - 1/Nguess*Tc;
    tb = tMax + 1/Nguess*Tc;
        
    hillaPrime = twoHillActivationDerivative( ta, m1, m2, tau1, tau2, Tc, Ts );
    hillbPrime = twoHillActivationDerivative( tb, m1, m2, tau1, tau2, Tc, Ts );
    
    % find critical point using bisection
    k = 0;
    while( max( abs(hillaPrime),abs(hillbPrime) ) > 1e-10 && k < 100)
        
        tc = (ta+tb)/2.0;
        hillcPrime = twoHillActivationDerivative( tc, m1, m2, tau1, tau2, Tc, Ts );
        
        
        if( hillaPrime >= 0 && hillcPrime <= 0 )
            tb = tc;
        else
            ta = tc;
        end
    
        hillaPrime = twoHillActivationDerivative( ta, m1, m2, tau1, tau2, Tc, Ts );
        hillbPrime = twoHillActivationDerivative( tb, m1, m2, tau1, tau2, Tc, Ts );
    
    
        k = k+1; 
    end
    
    tMax = (ta+tb)/2.0;
    hillMax = twoHillActivation( tMax, m1, m2, tau1, tau2, Tc, Ts, 1 );
   
 
   
end

