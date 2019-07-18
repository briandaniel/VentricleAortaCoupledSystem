function [ y_np1, U_LVOT, Z, w, U_AO ] = ODEstep( t, yn, U_LVOT, Z, w, U_AO, dt, param )
%ODEFUNC Summary of this function goes here
%   Detailed explanation goes here

    % Compute characteristics at new time step
    [wmInLVOT, wpOutLVOT] = computeOutgoingCharacteristics(U_LVOT, param.paramLVOT.x, dt, param.paramLVOT,2);
    [wmInAO, wpOutAO] = computeOutgoingCharacteristics(U_AO, param.paramAO.x, dt, param.paramAO,4);

        
    % Compute F at the previous time step
    [ Fn, Zn ] = Ffunc( yn, Z, wmInLVOT, t, param.paramLA, param.paramMV, param.paramLV, param );
    
    % Update the LA, MV, LV using a newton iteration.
    newtFunc = @(X,Z) Hfunc( X, Z, t, dt, Fn, yn, param, wmInLVOT );
    [ y_np1, Z ] = approximate_newton_iteration_Ndim_withZ( newtFunc, yn, Zn, param.newtParam );
    Plv = Z(4);

    % Compute lvot inflow boundary from previous information
    [ wpInLVOT ] = updateInflowLVOT( Plv, wmInLVOT, param );
   
    % Compute the flow through the aortic valve
    [ w, wmOutLVOT, wpInAO ] = updateAOV( w, dt, U_LVOT, U_AO, wpOutLVOT, wmInAO, param );
    
    % Compute the outflow from the AO
    Psa = w(3);
    [ Psa, wmOutAO ] = updateOutflowAO( Psa, dt, U_AO, wpOutAO, param);
    w(3) = Psa;
    
    % Update LVOT
    [ U1_LVOT, Uend_LVOT ] = computeBoundaries( wpInLVOT, wpOutLVOT, wmInLVOT, wmOutLVOT, param.paramLVOT );
    [ U_LVOT ] = vesselUpdateLW( U_LVOT, dt, param.paramLVOT.dx, param.paramLVOT, param.paramLVOT.x, param.paramLVOT.xmid, t );

    U_LVOT(1,:) = U1_LVOT;
    U_LVOT(end,:) = Uend_LVOT;
    
    
    % Update AO
    [ U1_AO, Uend_AO ] = computeBoundaries( wpInAO, wpOutAO, wmInAO, wmOutAO, param.paramAO );
    [ U_AO ] = vesselUpdateLW( U_AO, dt, param.paramAO.dx, param.paramAO, param.paramAO.x, param.paramAO.xmid, t );

    U_AO(1,:) = U1_AO;
    U_AO(end,:) = Uend_AO;
    
    
    
end

