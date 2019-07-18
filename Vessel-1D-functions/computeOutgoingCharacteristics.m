function [wmInlet, wpOutlet] = computeOutgoingCharacteristics(U, x, dt, paramVessel, m)
 
                                                            
    %Number of points to interpolate is m
    idx1 = 1:m;
    idx2 = (length(U(:,1))-m):length(U(:,1));
    
    % Inlet boundary
    lambdamInlet = U(idx1,2) - U(idx1,1).^(1/4).*paramVessel.R0(idx1).^(1/2);
    xInlet = x(idx1);
    R0Inlet = paramVessel.R0(idx1);
    WmInletPrev = U(idx1,2) - 4*sqrt(R0Inlet).*U(idx1,1).^(1/4) + 4*sqrt(R0Inlet).*paramVessel.A0(idx1).^(1/4); 
    wmInlet = interp1(xInlet + dt.*lambdamInlet,WmInletPrev,x(1));
        
    % Outlet boundary
    lambdapOutlet = U(idx2,2) + U(idx2,1).^(1/4).*paramVessel.R0(idx2).^(1/2);
    xOutlet = x(idx2);
    R0Outlet = paramVessel.R0(idx2);
    WpOutletPrev = U(idx2,2) + 4*sqrt(R0Outlet).*U(idx2,1).^(1/4) - 4*sqrt(R0Outlet).*paramVessel.A0(idx2).^(1/4); 
    wpOutlet = interp1(xOutlet + dt.*lambdapOutlet,WpOutletPrev,x(end));

end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

