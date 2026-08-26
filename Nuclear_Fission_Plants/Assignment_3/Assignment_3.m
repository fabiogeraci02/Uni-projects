clear; close all; clc;

%% Data

% Fuel pin geometry (AP1000)
inch_to_meter = 0.0254; % [m] to [in]
H_active = 168*inch_to_meter; % [m]  Active fuel height
D_co = 0.374*inch_to_meter; % [m]  Cladding outer-diameter
clad_thickness = 0.0225*inch_to_meter; % [m]  Cladding thickness
D_ci = D_co - 2*clad_thickness;% [m]  Cladding inner-diameter (cold)
zz = linspace(-H_active/2, H_active/2);

% Coolant condition
P_psia = 2250; % [psia]
P_bar = P_psia*0.0689476; % [bar]
P_sys = P_bar*1e5; % [Pa]

% Cladding temperature and pressure from assignmet 2
load('Tin_Tout.mat','Tin_Tout');
T_ci = max(Tin_Tout(1,:))+273.15; % [K]
T_ci_p = Tin_Tout(1, :); % [°C]
T_co_p = Tin_Tout(2, :); % [°C]
T_avg_p = (T_co_p+T_ci_p)/2+273.15; % [K]

load('P_T_in_cl.mat','P_T_in_cl');
p_cl_in = max(P_T_in_cl); % [Pa]

% Gas plenum
burn_up = 60000*1e6*86400/1e3; % [MW*Day/tons] to [J/kg]
fiss_yield = 0.28;
f_rate = 0.4;
Na = 6.02e23; % [-] Avogadro number
Rgas = 8.314; % [Pa*m^3/(mol*K)]

sigma_y = 241e6; % [Pa]
sigma_u = 413e6; % [Pa]

alpha_zr = @(TT) (6.72e-6*TT-2.07e-3)./(TT-308);
E_zry4 = @(TT) (9.9e3-5.669*(TT-273))*9.81e6; 
nu_cl = @(TT) 0.3303+8.376e-5*(TT-273);
k_cl_T = @(TT) 11.45+1.425e-2.*TT; 

%% 1) Preliminar Buckling

tt = clad_thickness;
E_ref = E_zry4(T_avg_p);
nu_ref = nu_cl(T_avg_p);

r_ci = D_ci/2;
r_co = D_co/2;
r_ave = 0.5*(r_ci+r_co);

% Critical pressure
p_cr = E_ref./(4*(1-nu_ref.^2))*(tt/r_ave)^3; % [Pa]
buck = p_cr/P_sys;

fprintf('\n1) PRELIMINARY BUCKLING CHECK\n');
fprintf('min p_cr  = %.2f  MPa\n', min(p_cr)/1e6);
fprintf('P_sys = %.2f  MPa  →  margin = %.1f ×\n', P_sys/1e6, min(buck));
if buck >= 1
    disp('Buckling margin OK (≥1)');
else
    warning('Insufficient buckling margin: reconsidering the design');
end

%% 2) Maximum internal pressure of the cladding with a depressurized reacto (p = 0) (sigma_h = sigma_y)

p_i_max = sigma_y*(tt/r_ave);

fprintf('2) Maximum internal pressure = %.2f MPa', p_i_max/1e6);

%% 3) Mechanical stresses

pin = p_i_max;
pout = P_sys;

Sm = min(2/3*sigma_y, 1/3*sigma_u);
% Inner wall stresses
s_h_in_m = (pin.*(r_co^2+r_ci^2)-2*pout*r_co^2)/(r_co^2-r_ci^2);
s_r_in_m = -pin;
s_l_in_m = r_ci^2.*(pin-pout)/(r_co^2-r_ci^2);

% Outer wall stresses
s_h_out_m = ((2*pin*r_ci^2)-pout.*(r_co^2+r_ci^2))/(r_co^2-r_ci^2);
s_r_out_m = -pout;
s_l_out_m = r_ci^2.*(pin-pout)/(r_co^2-r_ci^2);

% Inner Tresca criterium
s_max_in_m = max([abs(s_h_in_m - s_r_in_m),abs(s_l_in_m - s_r_in_m),abs(s_h_in_m - s_l_in_m)]);

% Outer Tresca criterium
s_max_out_m = max([abs(s_h_out_m - s_r_out_m),abs(s_l_out_m - s_r_out_m),abs(s_h_out_m - s_l_out_m)]);

% ASME verification
check_in  = s_max_in_m <= Sm;
check_out = s_max_out_m <= Sm;

fprintf('\n3) ASME CATEGORY P CHECK \n');
fprintf('S_m = %.2f MPa\n', Sm/1e6);
fprintf('Max σ_T (inner) = %.2f MPa → %s\n', s_max_in_m/1e6, tern(all(check_in),'OK','FAIL'));
fprintf('Max σ_T (outer) = %.2f MPa → %s\n', s_max_out_m/1e6, tern(all(check_out),'OK','FAIL'));

%% 4) Thermal stresses

T_ci_new = 417.80+273.15;
T_co_new = 347.28+273.15;
Tref_new = T_ci_new;
E_ref_new = E_zry4(Tref_new);
nu_ref_new = nu_cl(Tref_new);
alpha_ref_new = alpha_zr(Tref_new);


% Inner wall Thermal stresses
num_in = 1-2*r_co^2/(r_co^2-r_ci^2)*log(r_co/r_ci);
den_in = log(r_co/r_ci);

s_h_in_t = E_zry4(T_ci_new)*alpha_zr(T_ci_new)*(T_co_new-T_ci_new)/(2*(1-nu_cl(T_ci_new)))*num_in/den_in;
s_r_in_t = 0;
s_l_in_t = s_h_in_t;

% Outer wall thermal stresses
num_out = 1-2*r_ci^2/(r_co^2-r_ci^2)*log(r_co/r_ci);
den_out = log(r_co/r_ci);

s_h_out_t = E_zry4(T_co_new)*alpha_zr(T_co_new)*(T_co_new-T_ci_new)/(2*(1-nu_cl(T_co_new)))*num_out/den_out;
s_r_out_t = 0;
s_l_out_t = s_h_out_t;

% Inner thermal Tresca criterium
s_max_in_t = max([abs(s_h_in_t - s_r_in_t),abs(s_r_in_t - s_l_in_t),abs(s_h_in_t - s_l_in_t)]);

% Outer thermal Tresca criterium
s_max_out_t = max([abs(s_h_out_t - s_r_out_t),abs(s_r_out_t - s_l_out_t),abs(s_h_out_t - s_l_out_t)]);

% Primary + secondary inner
s_h_in = s_h_in_m+s_h_in_t;
s_r_in = s_r_in_m+s_r_in_t;
s_l_in = s_l_in_m+s_l_in_t;
% Primary + secondary outer
s_h_out = s_h_out_m+s_h_out_t;
s_r_out = s_r_out_m+s_r_out_t;
s_l_out = s_l_out_m+s_l_out_t;

% Tresca stresses P+Q 
tau_max_in = max([abs(s_h_in-s_r_in)/2,abs(s_r_in-s_l_in)/2,abs(s_h_in-s_l_in)/2]);
tau_max_out = max([abs(s_h_out-s_r_out)/2,abs(s_r_out-s_l_out)/2,abs(s_h_out-s_l_out)/2]);

%% ASME verification

limit_ASME = 3*Sm;
fprintf('\n5) [ASME CHECK] Limite: 3·Sm = %.2f MPa\n', limit_ASME/1e6);


% Inner wall
s_max_in = 2*tau_max_in;

fprintf('\n[INNER WALL]\nσ_max = %.2f MPa → ', s_max_in/1e6);
if s_max_in <= limit_ASME
    disp('OK: σ_max ≤ 3·Sm (ASME compliant)');
else
    disp('FAIL: σ_max > 3·Sm (not compliant)');
end

% Outer wall

s_max_out = 2*tau_max_out;

fprintf('\n[OUTER WALL]\nσ_max = %.2f MPa → ', s_max_out/1e6);
if s_max_out <= limit_ASME
    disp('OK: σ_max ≤ 3·Sm (ASME compliant)');
else
    disp('FAIL: σ_max > 3·Sm (not compliant)');
end


%% Size the gas plenum of the fuel pin

% Uranium dioxide mass
M_tot_UO2_lb = 211588;
M_tot_UO2_kg = M_tot_UO2_lb*0.45359237;
N_ass = 157;   
N_rod = 264; % rods of each assembly
N_tot = N_ass*N_rod; % total rods
M_UO2_pin = M_tot_UO2_kg/N_tot; % [kg] of each pin


m_N2 = 25e-6;
M_N2 = 28; % [g/mol]
m_H2O = 75e-6;
M_H2O = 18;

n_N2 = M_UO2_pin*m_N2/M_N2;
n_H2O = M_UO2_pin*m_H2O/M_H2O;


E_tot = burn_up*M_UO2_pin*(238/270);
E_fission = 200*1.6e-13;
N_fission = E_tot/E_fission; 

n_fiss = N_fission*fiss_yield*f_rate/Na; % [mol]
n_tot = n_fiss+n_H2O+n_N2;

V_pl = (n_tot*8.314*T_ci)/p_i_max;

H_pl = V_pl/(pi*r_ci^2)*100; % [cm]

fprintf('\n[Gas plenum size] = %.2f cm', H_pl); % without spring

% verification ideal gas hypotesis
P_h2o = (n_H2O/n_tot)*p_i_max;
Z_h2o = P_h2o*V_pl/(n_H2O*Rgas*Tref_new);

if Z_h2o>0.8
    fprintf('\nIdeal gas hypotesis verified');
else
    fprintf('\nIdeal gas hypotesis not verified');
end

fprintf('\nZ = %.3f \n', Z_h2o)

%% Post processing

% Buckling figure
figure(1)
yyaxis left
plot(zz, p_cr, 'LineWidth', 1.5)
ylabel('Critical pressure [Pa]')

yyaxis right
plot(zz, T_avg_p, 'LineWidth', 1.5)
ylabel('Temperature [K]')
grid on
xlabel('Core height [m]')
title('Critical pressure and Temperature profile in function of the position')
legend('Critical pressure', 'Average temperature')

% Table for inner stresses
data_in = [ s_h_in_m, s_r_in_m, s_l_in_m, s_max_in_m;
         s_h_in_t, s_r_in_t, s_l_in_t, s_max_in_t;
         s_h_in, s_r_in, s_l_in, 2*tau_max_in]/1e6;

rowLabels = {'Pm','Q','Pm+Q'};
colLabels = {'$\sigma_h$','$\sigma_r$','$\sigma_l$','$\sigma_{\mathrm{Tresca}}$'};


[nR, nC] = size(data_in);
totR = nR + 1;    
totC = nC + 1;    

figure('Color','w');
axis off
hold on

% vertical
for i = 0:totC
    if i==1
        lw = 2;     
    else
        lw = 1;
    end
    line([i i], [0 totR], 'Color','k', 'LineWidth', lw);
end
% orizzontal
for j = 0:totR
    if j==nR
        lw = 2;      
    else
        lw = 1;
    end
    line([0 totC], [j j], 'Color','k', 'LineWidth', lw);
end

for ic = 1:nC
    x = ic + 0.5;
    y = nR + 0.5;
    text(x, y, colLabels{ic},'Interpreter','latex', 'FontWeight','bold','HorizontalAlignment','center', 'VerticalAlignment','middle');
end

for jr = 1:nR
    x = 0.5;
    y = nR - jr + 0.5;
    text(x, y, rowLabels{jr},'FontWeight','bold','HorizontalAlignment','center', 'VerticalAlignment','middle');
end

for jr = 1:nR
    for ic = 1:nC
        x = ic + 0.5;
        y = nR - jr + 0.5;
        text(x, y, sprintf('%.2f MPa', data_in(jr,ic)),'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
end

xlim([0 totC]);
ylim([0 totR]);
axis equal
title('Inner stresses','FontWeight','bold');


% Table for outer stresses

data_out = [ s_h_out_m, s_r_out_m, s_l_out_m, s_max_out_m;
         s_h_out_t, s_r_out_t, s_l_out_t, s_max_out_t;
         s_h_out, s_r_out, s_l_out, 2*tau_max_out]/1e6;

rowLabels = {'Pm','Q','Pm+Q'};
colLabels = {'$\sigma_h$','$\sigma_r$','$\sigma_l$','$\sigma_{\mathrm{Tresca}}$'};


[nR, nC] = size(data_out);
totR = nR + 1;    
totC = nC + 1;    

figure('Color','w');
axis off
hold on

% Vertical
for i = 0:totC
    if i==1
        lw = 2;     
    else
        lw = 1;
    end
    line([i i], [0 totR], 'Color','k', 'LineWidth', lw);
end
% orizzontal
for j = 0:totR
    if j==nR
        lw = 2;      
    else
        lw = 1;
    end
    line([0 totC], [j j], 'Color','k', 'LineWidth', lw);
end

for ic = 1:nC
    x = ic + 0.5;
    y = nR + 0.5;
    text(x, y, colLabels{ic},'Interpreter','latex', 'FontWeight','bold','HorizontalAlignment','center', 'VerticalAlignment','middle');
end

for jr = 1:nR
    x = 0.5;
    y = nR - jr + 0.5;
    text(x, y, rowLabels{jr},'FontWeight','bold','HorizontalAlignment','center', 'VerticalAlignment','middle');
end

for jr = 1:nR
    for ic = 1:nC
        x = ic + 0.5;
        y = nR - jr + 0.5;
        text(x, y, sprintf('%.2f MPa', data_out(jr,ic)),'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
end

xlim([0 totC]);
ylim([0 totR]);
axis equal
title('Outer stresses','FontWeight','bold');

%% Function:

function out = tern(cond, val_true, val_false)
    if cond, out = val_true; else, out = val_false; end
end
