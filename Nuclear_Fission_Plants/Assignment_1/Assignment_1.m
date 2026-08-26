clear; close all; clc;

%% Data

P_nom = 6e8; % [W]
P_th = 0.01*P_nom;
dp_core = 1.2e5; % [Pa]
mcore = 3200; % [kg/s]

%PSC data
inch_to_meter = 0.0254;
p_c_PSC = 75e5; % [Pa]
D_out_PSC = 16*inch_to_meter; % [m]
L_PSC = 16; % [m]
eps_D_PSC = 2e-4;
kbend_PSC = 0.45;
Nbend_PSC = 4;
kloss_valve_PSC = 0.12;
H1_PSC = 7; % [m]
H2_PSC = 3; % [m]

D_in_PSC = (D_out_PSC-(2*1.031)*inch_to_meter);
A_pipe_PSC = pi*(D_in_PSC/2)^2;


%HX1 data
% Outer tube diameter [mm]
D_out_HX1 = 19.05e-3; % Converted in meter

% Tube thickness [mm]
t_HX1 = 1.24e-3; % Converted in meter

% Number of tubes
N_t_HX1 = 897;

A1_sec=pi*(D_out_HX1-2*t_HX1)^2/4;
A1_Nsec=A1_sec*N_t_HX1;

% Pitch [mm] 
d_p_HX1 = 28.5e-3; % Converted in meter

% Shell inner diameter [m]
D_shell_HX1 = 1.5;

% Number of baffles
N_b_HX1 = 2;

% Baffles spacing [m]
l_baffles_HX1 = 1.6;

% Average tube length [m]
L_tubes_HX1 = 9.314;

% Tube thermal conductivity [W/(m*K)]
k_th_HX1 = 15;

% Headers average flow area [m^2]
A_head_HX1 = 0.883;

% Heat transfer area [m^2]
A_HX1 = 500;

% Tube relative roughness [-]
eps_D_HX1 = 1e-4;

% Correction factor for ΔTml [-]
F_T = 0.7;

A1_in_sup = pi*(D_out_HX1-2*t_HX1)*L_tubes_HX1*N_t_HX1;
D_in_HX1=D_out_HX1-2*t_HX1;
%ISC data
p_c_ISC = 70e5; % [Pa]
D_out_ISC = 16*inch_to_meter; % [m]
L_ISC = 40; % [m]
eps_D_ISC = 2e-4;
kbend_ISC = 0.45;
Nbend_ISC = 6;
H_net_ISC = 10; % [m]
D_in_ISC = (D_out_ISC-(2*1.031)*inch_to_meter);
A_pipe_ISC = pi*(D_in_ISC/2)^2;

%HX2 data
% Tube outer diameter [mm] -> [m]
D_out_HX2 = 25.4e-3;

% Tube thickness [mm] -> [m]
t_HX2 = 1.24e-3;

% Number of tubes
N_t_HX2 = 770;

% Average tube length [m]
L_tubes_HX2 = 7;

% Manifold diameter [inc] -> [m]
D_manifold = 16 * 0.0254; % 1 inch = 0.0254 m

% Tube thermal conductivity [W/(m*K)]
k_th_HX2 = 15;

% Heat transfer area [m^2]
A_HX2 = 430;

% Tube relative roughness [-]
eps_D_HX2 = 1e-4;

%% Intermediate circuit

% Starting from HX2
q_HX2=P_th/A_HX2;
delta_sat=(q_HX2/2.257)^(1/3.86);
h_pool=q_HX2/delta_sat;

% Geometry HX2
D_in_HX2=D_out_HX2-2*t_HX2;
A2_sec=pi*(D_in_HX2)^2/4;
A2_Nsec = A2_sec*N_t_HX2;
A2_in_sup = pi*D_in_HX2*L_tubes_HX2*N_t_HX2;

% Geometry HX1
A_shell_HX1=(D_shell_HX1/d_p_HX1)*(d_p_HX1-D_out_HX1)*l_baffles_HX1;
Deq_shell=((2*sqrt(3)*d_p_HX1^2)/(pi*D_out_HX1))-D_out_HX1;

% Saturation temperature
T_sat = XSteam('tSat_p',1)+273.15; 

% Cicle starting point
Tm2 = 120+273.15; % [K]
m2 = 100; % [kg/s]

it2 = 0;
it2_m = 0;

itermax = 200;
itermax_m = 200;

err2 = 1000;
err2_m = 1000;

toll2 = 1e-6;
toll2_m = 1e-6;

g = 9.81;

while err2_m > toll2_m & it2_m < itermax_m
    it2_m = it2_m +1;
    
    it2=0;
    err2=1000;
    while err2 > toll2 & it2 < itermax
        it2 = it2+1; 
        
        % Properties of water
        rho2 = XSteam('rho_pt', 70, Tm2-273.15); % Density [kg/m³] a p = 70 bar e t = Tm2
        mu2 = XSteam('my_pt', 70, Tm2-273.15); % Dinamic viscosity [Pa·s]
        cp2 = XSteam('Cp_pt',70,Tm2-273.15)*1e3; % specific heat [J/(kg·K)]
        lambda2 = XSteam('tc_pt', 70,Tm2-273.15);  % Thermal conductivity [W/(m·K)]
        
        % Evaluation of the h_in
       
        vel2 = m2/(rho2*A2_Nsec);
        Re2 = rho2*vel2*D_in_HX2/mu2;
        fric2 = 1/(1.8*log10((eps_D_HX2/3.7)^1.11+6.9/Re2))^2;
    
        Pr2 = cp2*mu2/lambda2;
        
        Nu2 = ((fric2/8)*(Re2-1000)*Pr2)/(1+12.7*((fric2/8)^0.5)*(Pr2^(2/3)-1));
        h2 = lambda2*Nu2/D_in_HX2;
    
       
        U_HX2 = 1/(A_HX2*(1/(h2*A2_in_sup)+log(D_out_HX2/(D_out_HX2-2*t_HX2))/(2*pi*L_tubes_HX2*k_th_HX2*N_t_HX2)+1/(h_pool*A_HX2)));
    
        DeltaT_log2 = P_th/(U_HX2*A_HX2);
    
        aa=exp((U_HX2*A_HX2)/(m2*cp2));
        T_cold2=P_th/(m2*cp2*(aa-1))+T_sat;    
        T_hot2=T_cold2+P_th/(m2*cp2);
    
        Tm2_new = (T_hot2+T_cold2)/2;
    
        err2 = abs(Tm2_new-Tm2)/Tm2_new; 
    
        Tm2 = Tm2_new;
       
    end

    rho2_h = XSteam('rho_pt', 70, T_hot2-273.15);
    rho2_c = XSteam('rho_pt', 70, T_cold2-273.15);
    rho2 = XSteam('rho_pt', 70, Tm2-273.15);
    mu2 = XSteam('my_pt', 70, Tm2-273.15); % dinamic viscosity [Pa·s]
    
    vel_ISC = m2/(rho2*A_pipe_ISC);
    vel_shell=m2/(rho2*A_shell_HX1);
    
    Re_ISC = rho2*vel_ISC*(D_in_ISC)/mu2;
    Re_shell= rho2*vel_shell*Deq_shell/mu2;

    fric_ISC = 1/(1.8*log10((eps_D_ISC/3.7)^1.11+6.9/Re_ISC))^2;

    kshell=8*(0.227/(Re_shell^0.193))*(D_shell_HX1/Deq_shell)*((N_b_HX1+1)*l_baffles_HX1)/l_baffles_HX1;

    % Distributed losses
    dp_friction_HX2=fric2*(L_tubes_HX2/ D_in_HX2)/ (2 * rho2 * A2_Nsec^2);
    dp_friction_ISC=fric_ISC*(L_ISC/ D_in_ISC)/ (2 * rho2 * A_pipe_ISC^2);

    % Localized losses
    % Bend 
    dp_bend2=Nbend_ISC*kbend_ISC/(2*rho2*A_pipe_ISC^2);
    % Contraction hx2
    dp_contr_HX2=0.5*(1-(A2_sec/A_pipe_ISC))/(2*rho2_h*A2_Nsec^2);
    % Expansion hx2
    dp_exp_HX2=((1-(A2_sec/A_pipe_ISC))^2)/(2*rho2_c*A2_Nsec^2);
    % Contraction hx1
    dp_contr_HX1=0.5*(1-(A_pipe_ISC/A_shell_HX1))/(2*rho2_h*A_pipe_ISC^2);
    % Expansion hx1
    dp_exp_HX1=(1-(A_pipe_ISC/A_shell_HX1))^2/(2*rho2_c*A_pipe_ISC^2);
    
    % Shell
    dp_shell2=kshell/(2*rho2*A_shell_HX1^2);

    % Tot distributed losses
    tot_distr2=dp_friction_ISC+dp_friction_HX2;

    % Tot localized losses
    tot_local2=dp_bend2+dp_contr_HX2+dp_exp_HX2+dp_contr_HX1+dp_exp_HX1+dp_shell2;
   
    PI_2=tot_local2+tot_distr2;

    % Mass flow rate
    m2_new=sqrt(((rho2_c-rho2_h)*g*H_net_ISC)/PI_2);
    
    err2_m=abs(m2_new-m2)/m2_new;

    m2=m2_new;

end

%% Primary circuit

% Starting from HX1
rho2 = XSteam('rho_pt', 70, Tm2-273.15);
mu2 = XSteam('my_pt', 70, Tm2-273.15); 
cp2 = XSteam('Cp_pt',70,Tm2-273.15)*1e3; 
lambda2 = XSteam('tc_pt', 70,Tm2-273.15);  

Pr2 = cp2*mu2/lambda2;    
vel_shell=m2/(rho2*A_shell_HX1);
Re_shell= rho2*vel_shell*Deq_shell /mu2;
h_shell=0.351*Re_shell^0.55*(lambda2/Deq_shell)*Pr2^(1/3);

% Cycle starting point

m1 = 100; % [kg/s]

Tm1= 120+273.15; 

it1_m=0;
err1_m=1000;
toll1_m=1e-6;

while err1_m > toll1_m && it1_m < itermax_m
    it1_m = it1_m +1;
    
    it1=0;
    err1=1000;
    toll1=1e-6;
    while err1 > toll1 & it1 < itermax
            it1 = it1+1; 
        
        
        % Properties of water
        rho1 = XSteam('rho_pt', 75, Tm1-273.15); 
        mu1 = XSteam('my_pt', 75, Tm1-273.15); 
        cp1 = XSteam('Cp_pt',75, Tm1-273.15)*1e3; 
        lambda1 = XSteam('tc_pt', 75,Tm1-273.15); 
        
        % Evaluation of the h_in
        vel1 = m1/(rho1*A1_Nsec);
        Re1 = rho1*vel1*(D_in_HX1)/mu1;
        fric1 = 1/(1.8*log10((eps_D_HX1/3.7)^1.11+6.9/Re1))^2;
    
        Pr1 = cp1*mu1/lambda1;
        
        Nu1 = ((fric1/8)*(Re1-1000)*Pr1)/(1+12.7*((fric1/8)^0.5)*(Pr1^(2/3)-1));
        h1 = lambda1*Nu1/(D_in_HX1);
    
       
        U_HX1 = 1/(A_HX1*(1/(h1*A1_in_sup)+log(D_out_HX1/(D_in_HX1))/(2*pi*L_tubes_HX1*k_th_HX1*N_t_HX1)+1/(h_shell*A_HX1)));
    
        DeltaT_log1 = P_th/(U_HX1*A_HX1*F_T);
    
        bb=exp(U_HX1*A_HX1*F_T*(1/(m1*cp1) - (T_hot2-T_cold2)/P_th));
        T_cold1=(T_hot2-bb*T_cold2-P_th/(m1*cp1))/(1-bb);    
        T_hot1=T_cold1+P_th/(m1*cp1);
    
        Tm1_new = (T_hot1+T_cold1)/2;
    
        err1 = abs(Tm1_new-Tm1)/Tm1_new; 
    
        Tm1 = Tm1_new;
        
    end

     rho1_h = XSteam('rho_pt', 75, T_hot1-273.15);
     rho1_c = XSteam('rho_pt', 75, T_cold1-273.15);
     rho1 = XSteam('rho_pt', 75, Tm1-273.15);
     mu1 = XSteam('my_pt', 75, Tm1-273.15); % Dinamic viscosity [Pa·s]
     
     vel_PSC = m1/(rho1*A_pipe_PSC);
     Re_PSC = rho1*vel_PSC*D_in_PSC/mu1;
     fric_PSC = 1/(1.8*log10((eps_D_PSC/3.7)^1.11+6.9/Re_PSC))^2;
     
     % Distributed losses
     dp_friction_HX1=fric1*(L_tubes_HX1/ D_in_HX1)/ (2 * rho1 * A1_Nsec^2);
     dp_friction_PSC=fric_PSC*(L_PSC/ D_in_PSC)/ (2 * rho1 * A_pipe_PSC^2);
     
     % Tot distributed losses
     tot_distr1 = dp_friction_PSC+dp_friction_HX1;
     
     % Localized losses
     % Contraction header-U pipe
     dp_contr_header_HX1=0.5*(1-(A1_sec/A_head_HX1))/(2*rho1_h*A1_Nsec^2);
     % Expansion U pipe-header
     dp_exp_HX1_header=((1-(A1_sec/A_head_HX1))^2)/(2*rho1_c*A1_Nsec^2);

     % Contraction header-PSC 
     dp_contr_header_PSC=0.5*(1-(A_pipe_PSC/A_head_HX1))/(2*rho1_c*A_pipe_PSC^2);
     % Expansion PSC-header
     dp_exp_PSC_header=(1-(A_pipe_PSC/A_head_HX1))^2/(2*rho1_h*A_pipe_PSC^2);

     % Bend
     dp_bend1=Nbend_PSC*kbend_PSC/(2*rho1*A_pipe_PSC^2);

     % Valve
     dp_valve=kloss_valve_PSC/(2*rho1_c*A_pipe_PSC^2);
     
     % Vessel
     dp_vessel=dp_core/mcore^2;

     % Tot localized losses
     tot_local1=dp_bend1+dp_contr_header_HX1+dp_exp_HX1_header+dp_contr_header_PSC+dp_exp_PSC_header+dp_valve+dp_vessel;
   
     PI_1=tot_local1+tot_distr1;

     % Buoyancy forces
     dp_core_buoy = -g*(rho1_c+rho1_h)*H2_PSC/2;
     dp_PSC_buoy_c = g*rho1_c*(H1_PSC+H2_PSC);
     dp_PSC_buoy_h = -g*rho1_h*H1_PSC;
    
     
     % Flow
     m1_new=sqrt((dp_PSC_buoy_h + dp_PSC_buoy_c + dp_core_buoy)/PI_1);

     err1_m = abs(m1_new-m1)/m1_new;

     m1 = m1_new;
end


T_cold2 = T_cold2-273.15;
T_hot2 = T_hot2-273.15;
T_cold1 = T_cold1-273.15;
T_hot1 = T_hot1-273.15;

%% Table for ISC and PSC

% ISC
dataISC = {
    'Hot Leg Temperature (HX2)', sprintf('%.2f °C', T_hot2);
    'Cold Leg Temperature (HX2)', sprintf('%.2f °C', T_cold2);
    'Mass Flow Rate (ISC)', sprintf('%.2f kg/s', m2);
    'Distributed loss ISC', sprintf('%.2f Pa', dp_friction_ISC * m2^2);
    'Distributed loss HX2', sprintf('%.2f Pa', dp_friction_HX2 * m2^2);
    'Bend Loss', sprintf('%.2f Pa', dp_bend2 * m2^2);
    'Contraction HX2', sprintf('%.2f Pa', dp_contr_HX2 * m2^2);
    'Expansion HX2', sprintf('%.2f Pa', dp_exp_HX2 * m2^2);
    'Contraction HX1', sprintf('%.2f Pa', dp_contr_HX1 * m2^2);
    'Expansion HX1', sprintf('%.2f Pa', dp_exp_HX1 * m2^2);
    'Shell Loss', sprintf('%.2f Pa', dp_shell2 * m2^2);
    'Total Localized Loss', sprintf('%.2f Pa', tot_local2 * m2^2);
    'Total Distributed Loss', sprintf('%.2f Pa', tot_distr2 * m2^2);
    'Total Pressure Drop ISC', sprintf('%.2f Pa', PI_2 * m2^2)
};

fig1 = figure('Name', 'ISC Summary', 'NumberTitle', 'off', 'Color', [1 1 1]);
uicontrol('Style', 'text','String', 'ISC Summary','FontSize', 14,'FontWeight', 'bold','BackgroundColor', [1 1 1],'Position', [20 405 480 25]);

uitable('Data', dataISC,'ColumnName', {'Description', 'Value'},'ColumnWidth', {300, 150},'FontSize', 12,'RowName', [],'Position', [20 20 480 380]);

% PSC
dataPSC = {
    'Hot Leg Temperature (HX1)', sprintf('%.2f °C', T_hot1);
    'Cold Leg Temperature (HX1)', sprintf('%.2f °C', T_cold1);
    'Mass Flow Rate (PSC)', sprintf('%.2f kg/s', m1);
    'Distributed loss PSC', sprintf('%.2f Pa', dp_friction_PSC * m1^2);
    'Distributed loss HX1', sprintf('%.2f Pa', dp_friction_HX1 * m1^2);
    'Contraction header-HX1', sprintf('%.2f Pa', dp_contr_header_HX1 * m1^2);
    'Expansion HX1-header', sprintf('%.2f Pa', dp_exp_HX1_header * m1^2);
    'Contraction header-PSC', sprintf('%.2f Pa', dp_contr_header_PSC * m1^2);
    'Expansion PSC-header', sprintf('%.2f Pa', dp_exp_PSC_header * m1^2);
    'Bend Loss', sprintf('%.2f Pa', dp_bend1 * m1^2);
    'Valve Loss', sprintf('%.2f Pa', dp_valve * m1^2);
    'Vessel Loss', sprintf('%.2f Pa', dp_vessel * m1^2);
    'Total Localized Loss', sprintf('%.2f Pa', tot_local1 * m1^2);
    'Total Distributed Loss', sprintf('%.2f Pa', tot_distr1 * m1^2);
    'Total Pressure Drop PSC', sprintf('%.2f Pa', PI_1 * m1^2)
};

fig2 = figure('Name', 'PSC Summary', 'NumberTitle', 'off', 'Color', [1 1 1]);
uicontrol('Style', 'text','String', 'PSC Summary','FontSize', 14,'FontWeight', 'bold', 'BackgroundColor', [1 1 1],'Position', [20 405 480 25]);
uitable('Data', dataPSC,'ColumnName', {'Description', 'Value'},'ColumnWidth', {300, 150},'FontSize', 12,'RowName', [],'Position', [20 20 480 380]);

exportgraphics(fig1, 'ISC_Table.pdf', 'ContentType', 'vector');
exportgraphics(fig2, 'PSC_Table.pdf', 'ContentType', 'vector');