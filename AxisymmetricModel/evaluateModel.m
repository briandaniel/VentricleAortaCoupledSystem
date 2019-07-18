function [ alpha, kappa, eta, dVda, Vlv ] = evaluateModel( a1, a2, a3, t, paramLV )
%EVALUATEMODEL Summary of this function goes here
%   Detailed explanation goes here

    % Calculate new mu values
    muvec0 = linspace(paramLV.muin0,paramLV.muout0,paramLV.Nmu);
    nu_vec = linspace(paramLV.nu_up,pi,paramLV.Nnu);

    mu = zeros(paramLV.Nmu,paramLV.Nnu);
    for k=1:paramLV.Nmu
        mu(k,:) = muCalc( paramLV.a0, a1, a2, muvec0(k), nu_vec, paramLV.muin0);
    end

    subel = subsub_elements( mu, nu_vec, paramLV.a0, a1, a2, muvec0, paramLV.muin0 );
    subel.muvec0 = muvec0;

    % Test deformation tensor
    [ F, dF_da1, dF_da2, dF_da3 ] ...
        = deformation_tensor_derivatives( subel, paramLV, a1, a2, a3);

    % Test cauchy green defrormation tensor
    [ C ] = cauchy_tensor( F );

    % Test green strain tensor
    [ E, dE_da1, dE_da2, dE_da3 ] ...
                 = green_strain_derivatives( F, dF_da1, dF_da2, dF_da3, C );

    % Generate rotation matrix
    [ cps, sps, omega ] = initial_rotation_matrix( subel, paramLV );

    % generate determinant of F
    [ detF ] = determinantF( F );

    % Generate Cinv
    [ Cinv ] = inverse_cauchy_green_deformation( C );

    % Generate Sv the viscous stress
    [ Sv_a1 ] = viscous_stress( paramLV.kv, detF, Cinv, dE_da1 );
    [ Sv_a2 ] = viscous_stress( paramLV.kv, detF, Cinv, dE_da2 );
    [ Sv_a3 ] = viscous_stress( paramLV.kv, detF, Cinv, dE_da3 );

    % Generate Se the elastic stress
    [ Se ] = elastic_stress ( C, E, sps, cps, paramLV );

    % Generate derivatives of ell
    [ ell, ell0, dell_da1, dell_da2, dell_da3 ]  = ...
                ell_derivative( subel, paramLV, paramLV.a0, a1, a2, a3, omega );

    % Compute end-diastolic fiber strain 
    if(paramLV.useTwoHill == 1)
        At = twoHillActivation( t, paramLV.m1, paramLV.m2, paramLV.tau1, paramLV.tau2, paramLV.Tc, paramLV.Ts, paramLV.hillMaxVal );
    else
        At = activation_func( t, paramLV.Ta, paramLV.Tc, paramLV.eps_fed, paramLV.kp );
    end
    
    
    % Active fiber stress function  
    [ Sf_const ] = active_fiber_stress_constant( At, paramLV.km, paramLV.ked, paramLV.m, ...
                                         ell, ell0, paramLV.eps_fed, sps, cps, F );
    [ Sf_a1 ] = active_fiber_stress( At, paramLV.kav, paramLV.ked, paramLV.m, ...
                                         ell, ell0, paramLV.eps_fed, dell_da1, sps, cps, F );
    [ Sf_a2 ] = active_fiber_stress( At, paramLV.kav, paramLV.ked, paramLV.m, ...
                                         ell, ell0, paramLV.eps_fed, dell_da2, sps, cps, F );
    [ Sf_a3 ] = active_fiber_stress( At, paramLV.kav, paramLV.ked, paramLV.m, ...
                                         ell, ell0, paramLV.eps_fed, dell_da3, sps, cps, F );

    % Traction integral
    [eta1, eta2] = traction_integral( subel, paramLV, nu_vec, paramLV.a0, a1, a2  );

    % Left hand equilibrium integral
    [ alpha_a1, beta_a1, gamma_a1, kappa_a1 ] ...
        = equilibrium_integral( paramLV.a0, a1, muvec0, nu_vec, dE_da1, Sv_a1, Sv_a2, Sv_a3, ...
                                 Sf_const, Sf_a1, Sf_a2, Sf_a3, Se);
    [ alpha_a2, beta_a2, gamma_a2, kappa_a2 ] ...
        = equilibrium_integral( paramLV.a0, a1, muvec0, nu_vec, dE_da2, Sv_a1, Sv_a2, Sv_a3, ...
        Sf_const, Sf_a1, Sf_a2, Sf_a3, Se);
    [ alpha_a3, beta_a3, gamma_a3, kappa_a3 ] ...
        = equilibrium_integral( paramLV.a0, a1, muvec0, nu_vec, dE_da3, Sv_a1, Sv_a2, Sv_a3, ...
        Sf_const, Sf_a1, Sf_a2, Sf_a3, Se);

    % Calculate volume derivatives
    [ dV_da1, dV_da2 ] = volume_derivatives( subel, paramLV);
   
    % Calculate LV volume
    [ Vlv ] = cavity_volume_calc(subel.mu_in, nu_vec, subel.dmu_dnu_in, subel.a);
 
    % Turn these values into vectors/matrices as defined
    alpha = [alpha_a1, beta_a1, gamma_a1; ...
             alpha_a2, beta_a2, gamma_a2;
             alpha_a3, beta_a3, gamma_a3];
    kappa = [kappa_a1; kappa_a2; kappa_a3];
    eta = [eta1; eta2; 0];
    dVda = [dV_da1; dV_da2; 0];
    

end











