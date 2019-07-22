
clear all
close all


% run modelParameters1.m
% run LVparameters1.m

run modelParameters2.m
run LVparameters2.m

% run modelParameters3.m
% run LVparameters3.m


%%

imgPath = 'Images/';
imgIdx = 1;

% Time domain
Nc = 10;
tmax = param.T*Nc;
dt = 0.0005;
t = 0;

Nstep = tmax/dt;
plotStepSkip = 10;
textStepSkip = 10;

% Initial values
a1 = 0;
a2 = 0;
a3 = 0;
Pla = 0;
qmv = 0; 
zetamv = 0;
qaov = 0;
zetaaov = 0;
Psa = 80*1333;


% LVOT initial conditions
u0_LVOT = zeros(paramLVOT.Nx,1);
A0_LVOT = paramLVOT.A0;
U0_LVOT = [ A0_LVOT, u0_LVOT ];
U_LVOT = U0_LVOT;

% AO initial conditions
u0_AO = zeros(paramAO.Nx,1);
P0_AO = Psa*ones(paramAO.Nx,1);
A0_AO = paramAO.A0.*(P0_AO./paramAO.alpha + 1).^2;
U0_AO = [ A0_AO, u0_AO ];
U_AO = U0_AO;
P_AO = paramAO.alpha.*( sqrt(U_AO(:,1)./paramAO.A0) - 1 );

% These should probably approximated before solution loop
Z = [0, 0, 0, 0]'; % Starting guess for derivatives, etc.
Plv = Z(4);
y0 = [ a1, a2, a3, Pla, qmv, zetamv ]';
y = y0;
w = [qaov, zetaaov, Psa]';




%%
Vlv0 = cavity_volume_calc_minimum_input( a1, a2, paramLV );
Vlv = Vlv0;

% History
tvec = t;
Plv_vec = Plv;
Pla_vec = Pla;
qmv_vec = qmv;
zetamv_vec = zetamv;
qaov_vec = qaov;
zetaaov_vec = zetaaov;
a1_vec = a1;
a2_vec = a2;
a3_vec = a3;
Vlv_vec = Vlv;
Psa_vec = Psa;
Pao_end_vec = P_AO(end);
Pao_start_vec = P_AO(1);

% Check variables
Vla_approx = Pla*paramLA.C;
Ventered = 0.0;
Vexited = 0.0;
Vqmv = 0.0;
VinLa = 0.0;

Vla0 = Pla*paramLA.C;
Vsa0 = Psa*paramSA.C;
Vlvot0 = trapz( paramLVOT.x, U_LVOT(:,1) );
Vao0 = trapz( paramAO.x, U_AO(:,1) );
Vsystem0 = Vla0 + Vlv0 + Vlvot0 + Vao0 + Vsa0;

Vlvot_vec = Vlvot0;

%%
close all
figure('outerposition',[ 500, 25, 1600, 1000 ])
sp1a  = subplot('position',[0.05,0.8,0.25,0.15]);
sp1b  = subplot('position',[0.05,0.6,0.25,0.15]);
sp2 = subplot('position',[0.35,0.6,0.6,0.35]);
sp3a = subplot('position',[0.05,0.3,0.25,0.18]);
sp3b = subplot('position',[0.05,0.1,0.25,0.18]);
sp4 = subplot('position',[0.37,0.1,0.25,0.4]);
sp5 = subplot('position',[0.69,0.1,0.25,0.4]);


%%
tic
for k = 2:Nstep

    [ y, U_LVOT, Z, w, U_AO ] = ODEstep( t, y, U_LVOT, Z, w, U_AO, dt, param );

    t = t+dt;
        
    Plv = Z(4);
    a1 = y(1);
    a2 = y(2);
    a3 = y(3);
    Pla = y(4);
    qmv = y(5);
    zetamv = y(6);
    qaov = w(1);
    zetaaov = w(2);
    Psa = w(3);
    Vlv = cavity_volume_calc_minimum_input( a1, a2, paramLV );
    P_AO = paramAO.alpha.*( sqrt(U_AO(:,1)./paramAO.A0) - 1 );

    
    a1_vec = [a1_vec;a1];
    a2_vec = [a2_vec;a2];
    a3_vec = [a3_vec;a3];
    Plv_vec = [ Plv_vec; Plv ];
    Pla_vec = [ Pla_vec; Pla ];
    Psa_vec = [ Psa_vec; Psa ];
    qmv_vec = [ qmv_vec; qmv ];
    zetamv_vec = [ zetamv_vec; zetamv ];
    Vlv_vec = [Vlv_vec;Vlv];
    qaov_vec = [ qaov_vec; qaov ];
    zetaaov_vec = [ zetaaov_vec; zetaaov ];
    Pao_end_vec = [ Pao_end_vec ; P_AO(end) ];
    Pao_start_vec = [ Pao_start_vec; P_AO(1) ];
    tvec = [ tvec; t ];
  
    Vlvot = trapz( paramLVOT.x, U_LVOT(:,1) );
    Vlvot_vec = [Vlvot_vec;Vlvot];

    % Check computations
    % Inflow/outflow
    qin_la = (param.Ppv - Pla)/paramLA.R;
    qout_sa = (Psa - param.Psp)/paramSA.R;
        
    Ventered = Ventered + qin_la*dt;
    Vexited = Vexited + qout_sa*dt; 
    
    Vla = Pla*paramLA.C;
    Vsa = Psa*paramSA.C;
    Vao = trapz( paramAO.x, U_AO(:,1) );
    Vsystem = Vla + Vlv + Vlvot + Vao + Vsa;
        
    if( mod(k,textStepSkip) == 0 )
        display(['t = ', num2str(t), ',   Venter - Vexit = ', num2str(Ventered-Vexited), ...
            ',   Vsystem = ', num2str(Vsystem ), ',    VsystemDiff = ', num2str(Vsystem - Vsystem0)] )
    end
    
    % Plot
    if( mod(k,plotStepSkip) == 0 )
        
        subplot(sp1a)
            plot(tvec, qmv_vec, tvec, qaov_vec)
            ylabel('q_valve [mL/s]')
            title(['k = ',num2str(k),'    t = ', num2str(t)])
            xlabel('t [s]')
            
        subplot(sp1b)
            plot(tvec, zetamv_vec, tvec, zetaaov_vec)
            ylabel('\zeta_valve')
            xlabel('t [s]')
            
        subplot(sp2)
%             plot(param.paramLVOT.x,U_LVOT(:,2),...
%                  param.paramAO.x+param.paramLVOT.x(end),U_AO(:,2))
%             title(['k = ',num2str(k),'    t = ', num2str(t)])
%             ylabel('u [cm/s]')
%             Plvot = param.paramLVOT.alpha.*( sqrt( U_LVOT(:,1)./ param.paramLVOT.A0 ) - 1 );
%             Pao = param.paramAO.alpha.*( sqrt( U_AO(:,1)./ param.paramAO.A0 ) - 1 );
%             plot(param.paramLVOT.x,Plvot/1333,param.paramAO.x+param.paramLVOT.x(end),Pao/1333)
plot(param.paramAO.x,param.paramAO.A0,param.paramAO.x,U_AO(:,1))
            xlabel('x [cm]')
            set(gca,'ylim',[-10,180])
            set(gca,'xlim',[0,param.paramLVOT.x(end)+param.paramAO.x(end)])
            ylabel('P [mmHg]')

        subplot(sp3a)
            plot( tvec, Plv_vec/1333, tvec, Pla_vec/1333, tvec, Psa_vec/1333)
            ylabel('P [mmHg]')
            
        subplot(sp3b)
            plot(tvec, Pao_start_vec/1333, tvec, Pao_end_vec/1333, tvec, Psa_vec/1333);
            ylabel('P [mmHg]')
            xlabel('t [s]')

        subplot(sp4)
%             wp_ao = U_AO(:,2) + 4*paramAO.R0.^(1/2).*U_AO(:,1).^(1/4) - 4.*paramAO.R0.^(1/2).*paramAO.A0.^(1/4);
%             wm_ao = U_AO(:,2) - 4*paramAO.R0.^(1/2).*U_AO(:,1).^(1/4) + 4.*paramAO.R0.^(1/2).*paramAO.A0.^(1/4);
%             plot(param.paramAO.x, wp_ao, param.paramAO.x, wm_ao)
%             ylabel('V [mL]')
%             xlabel('t [s]')
        plot(Vlv_vec + Vlvot_vec,Plv_vec/1333)
        
        subplot(sp5)
            plot( tvec, a1_vec, tvec, a2_vec, tvec, a3_vec );
            ylabel('a_i')
            xlabel('t [s]')
            
            
        img = getframe(gcf);
%         display('Writing image to file...')
%         display([imgPath,'img',num2str(imgIdx),'.png'])
%         imwrite(img.cdata,[imgPath,'img',num2str(imgIdx),'.png']);
%         imgIdx = imgIdx+1;


    end

end
toc




save('coupled_simulation.mat')




































