function [ x,Z ] = approximate_newton_iteration_Ndim_withZ( ffunc, x0, Z0, newtParam )
%APPROXIMATE_NEWTON_ITERATION Summary of this function goes here
%   Detailed explanation goes here
    
    dx = newtParam.dx;
    xMinDiff = newtParam.xMinDiff;
    fNewtMin = newtParam.fNewtMin;
    maxIter = newtParam.maxIter;


    N = length(x0);

    x = x0;
    Z = Z0;
    xdiff = xMinDiff*2;
    k = 1;
    fx = fNewtMin*2;
    while ( norm(xdiff) > xMinDiff) && (norm(fx) > fNewtMin && (k <= maxIter) )
        
        [fx,Z] = ffunc(x, Z);
        
%         plot3(x(1),x(2),fx(1),'r.','markersize',16)
%         
%        display(['at x = ',num2str(x'),', fx = ',num2str(fx')])
        [ J ] = approximateJacobian( ffunc, fx, x, N, dx, Z );
        
        xdiff = J\fx;
        
        x = x - xdiff;

        k = k+1;
        
    end


end



function [ J ] = approximateJacobian( f, fx, x, N, dx, Z )

    J = zeros(N,N);

    for k = 1:N
        xdxk = x;
        
        xdxk(k) = x(k) + dx;
        fdxk = f(xdxk,Z);
        J(:,k) = ( fdxk - fx )./dx; % Compute approximat derivative
       
    end

end