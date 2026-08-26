clear; close all; clc;

%% Data
inch_to_meter = 0.0254; % [inc] to [m]

P_th = 3400e6; % [W]
P_th_core=P_th*0.974; % [W]
H_active = 168*inch_to_meter; % [m]
N_f_ass = 157; % [-]
N_f_rod = 264; % [-]
Ntot_f_rod = N_f_rod*N_f_ass; % [-]

D_co = 0.374*inch_to_meter; % [m]
clad_thickness = 0.0225*inch_to_meter; % [m]
D_ci = D_co-2*clad_thickness; % [m]
Dgap = 0.0065*inch_to_meter; % [m]
D_p = 0.3225*inch_to_meter; % [m]
A_sec_f = pi*D_p^2/4; % [m^2]
V_f_rod = A_sec_f*H_active; % [m^3]
Vtot_f_rod = Ntot_f_rod*V_f_rod; % [m^3]

% Point 2
F_Q = 2.6; % [-]
lambda_tr = 0.29e-2; % [m]
delta_ref = 1.72e-2; % [m]

H_e = H_active+1.42*lambda_tr+2*delta_ref; % [m]

% Point 3
conv_mass = 0.4535924/3600; % [lb/h] to [kg/s]
conv_area = 0.092903; % [ft^2] to [m^2]
mflow_eff = 106.8e6*conv_mass; % [kg/s]
area_eff = 41.8*conv_area; % [m^2]

% Point 4
T_in_F = 535; % [F]
T_in_celsius = (T_in_F-32)*5/9; % [°C]
T_in =T_in_celsius+273.15; % [K]
P_psia = 2250; % [psia]
P_bar = P_psia * 0.0689476; % [bar]
r_pitch = 0.496*inch_to_meter; % [m]

%% 1 Average volumetric heat generation rate

q_vol_av = P_th_core/Vtot_f_rod;

fprintf('Average volumetric heat generation rate = %.2f MW/m^3\n', q_vol_av*1e-6);

%% 2 Maximum volumetric heat generation rate

q_vol_max = q_vol_av*F_Q;

fprintf('Maximum volumetric heat generation rate = %.2f MW/m^3\n', q_vol_max*1e-6);

%% 3 Average mass velocity

G_ave = mflow_eff/area_eff;

fprintf('Average mass velocity = %.2e kg/(s*m^2)\n', G_ave);

%% 4 Coolant specific enthalpy and temperature profile

i_in = XSteam('h_pT', P_bar, T_in_celsius);
T_sat = XSteam('TSat_p',P_bar);
T_ave = (T_sat + T_in_celsius)/2;
cp_ave = XSteam('Cp_pT',P_bar,T_ave)*1e3;

A_subch = r_pitch^2 - pi*D_co^2/4;
mflow_subch = G_ave*A_subch;
zz = linspace(-H_active/2,H_active/2);
i_z = i_in + (1.0267*q_vol_max*A_sec_f*H_e/(mflow_subch*pi)*(sin(pi*zz/H_e)+sin(pi*H_active/(2*H_e))))*1e-3;

figure(1)
plot(i_z, zz,LineWidth=2)
grid on
xlabel("Specific enthalpy [kJ/kg]")
ylabel('Core height [m]')
title("Axial specific enthalpy profile")

% Tw = T_in_celsius + (1.0267*q_vol_max*A_sec_f*H_e/(mflow_subch*pi*cp_ave)*(sin(pi*zz/H_e)+sin(pi*H_active/(2*H_e))));
% Tw(Tw >= T_sat) = T_sat;
Tw = XSteamW('T_ph',P_bar,i_z);
figure(2)
plot(Tw,zz,LineWidth=2)
hold on 
xline(T_sat,'--','Saturation temperature','LineWidth',1.5,'LabelVerticalAlignment','bottom')
grid on
xlabel("Temperature [°C]")
ylabel('Core height [m]')
title("Coolant temperature profile")

%% 5 Equilibrium quality profile

i_lsat = XSteam('hL_p',P_bar);
i_vsat = XSteam('hV_p',P_bar);

x_eq = (i_z - i_lsat)/(i_vsat-i_lsat);

figure(3)
plot(x_eq,zz,LineWidth=1.5)
grid on
xlabel('Quality [kgv/kg]')
ylabel('Core height [m]')
title('Equilibrium quality profile')

%% 6.1 Cladding outer-wall temperature profile assuming negligible the presence of partial subcooled boiling region

% Convective heat transfer coefficient in sigle-phase
pH = pi*D_co; % heated perimeter
D_eq = 4*A_subch/pH; % equivalent diameter of the subchannel
CC = 0.042*r_pitch/D_co-0.024; % coefficient for the subchannel geometry
    
mu_w = zeros(length(Tw),1);
cp_w = zeros(length(Tw),1);
k_w = zeros(length(Tw),1);
Re = zeros(length(Tw),1);
Pr = zeros(length(Tw),1);
Nu = zeros(length(Tw),1);
hL = zeros(length(Tw),1);

for ii = 1:length(Tw)
    if Tw(ii) < T_sat

        mu_w(ii) = XSteam('my_pT', P_bar,Tw(ii));
        cp_w(ii) = XSteam('cp_ph', P_bar,Tw(ii))*1e3;
        k_w(ii) = XSteam('tc_pT', P_bar,Tw(ii));
    else
        mu_w(ii) = mu_w(ii-1);
        cp_w(ii) = cp_w(ii-1);
        k_w(ii) = k_w(ii-1);
    end    
   
    % Dittus-Boelter equation for turbulent flow
    Re(ii) = G_ave*D_eq/mu_w(ii);
    Pr(ii) = cp_w(ii)*mu_w(ii)/k_w(ii);
    Nu(ii) = CC*Re(ii)^0.8*Pr(ii)^0.4;
    hL(ii) = Nu(ii)*k_w(ii)/D_eq;
end

% T_co without subcooled boiling 
q_vol = q_vol_max*cos(pi*zz/H_e);
q_flux = q_vol.*D_p^2./(4*D_co);
T_co_SP = Tw + q_flux./hL';

% Jens-Lottes correlation
T_co_JL = T_sat + 25*(q_flux*1e-6).^(0.25)*exp(-P_bar/62);

% Griffith model
T_D_Gr = T_sat - q_flux./(5*hL');

% Position of ONB and detachment
zONB = 0;
zD = 0;
for yy = 1:length(Tw) 
    if T_co_JL(yy) <= T_co_SP(yy) && zONB == 0
        zFDB = zz(yy:end);
        zONB = zFDB(1);
        yONB = yy;
        fprintf("ONB a z = %.2f m\n", zONB+H_active/2)
    end

    if Tw(yy) >= T_D_Gr(yy) && zD == 0
        zD = zz(yy); % beginning of the detachment
        yD = yy;
        fprintf('Detachment at position: %.2f m\n', zD+H_active/2)
    end
end

% T_co considering subcooled boiling
T_co_sub_JL = T_co_SP;
T_co_sub_JL(zz >= zONB) = T_co_JL(zz >= zONB);

% plot with J-L correlation
figure(4)
plot(T_co_SP,zz,LineWidth=1.5)
hold on
plot(T_co_JL,zz,LineWidth=1.5)
plot(T_co_sub_JL,zz,LineWidth=1.5,LineStyle="--",Color='g')
grid on
xlabel('Temperature [°C]')
ylabel('Core height [m]')
title('Cladding outer surface temperature profile')
legend('T_{S-P}','T_{J-L}','T_{co}',Location='northwest')

figure(5)
plot(T_co_sub_JL,zz,LineWidth=1.5,Color='r')
hold on
plot(Tw,zz,LineWidth=1.5,Color='b')
grid on
xlabel('Temperature [°C]')
ylabel('Core height [m]')
legend('outer-cladding','coolant',Location='northwest')

%% 6.2 Flow quality profile

% Bowring-Rouhani model
rhoL_w = XSteamW('rho_pT', P_bar, Tw(Tw<T_sat));
rho_w_vsat = XSteam('rhoV_T', T_sat);
rho_w_lsat = XSteam('rhoL_T',T_sat);
rho_wsat = 1 ./ ((1 - x_eq(Tw == T_sat)) / rho_w_lsat + x_eq(Tw == T_sat) / rho_w_vsat);
rho_wtot = [rhoL_w, rho_wsat];
i_evap = i_vsat-i_lsat; % [kJ/kg] evaporation enthalpy
i_lp = i_z;
i_lp(i_z >= i_lsat) = (i_z(i_z >= i_lsat) - i_vsat*x_eq(i_z >= i_lsat))./(1-x_eq(i_z >= i_lsat)); % liquid phase enthalpy
eps = rho_wtot.*(i_lsat-i_lp)/(rho_w_vsat*i_evap);
eps(eps < 0) = 0;
zD_out = zz(yD:end);
xx = pH*q_vol_max*D_p^2./(i_evap*1e3*mflow_subch.*(1+eps(zz >= zD))*4*D_co).*(sin(pi*zD_out/H_e) - sin(pi*zD/H_e));

figure(6)
plot(xx,zD_out + H_active/2,'LineWidth',1.5)
hold on
plot(x_eq(yD:end),zD_out + H_active/2,'LineWidth',1.5)
grid on
title('Flow quality profile')
xlabel('Flow quality [kgv/kg]')
ylabel('Core height [m]')
legend('B-R model','Equilibrium quality','Location','northwest')

%% 6.3 Void fraction profile

Rb = 2.37e-3/P_bar^(0.237); % bubble radius by Rouhani
delta_b = 0.0666*Rb; % bubble layer tickness by Bowring

% Void fraction at the detachment by Maurer
alpha_D = 4*delta_b/D_eq;

% Void fraction profile from detachment to outlet by Zuber-Findlay correlation
% Constant by Collier-Thome
C0 = 1.13;
C1 = 1.18;

Tc = 647.096;  % critical temperature of water (K)
Tw_K = Tw + 273.15;  
sigma = 0.2358 * (1 - Tw_K/Tc).^1.256;  % surface tension in N/m

alphaD_out = xx/rho_w_vsat./(C0*(xx/rho_w_vsat + (1-xx)/rho_w_lsat) + C1./G_ave*(sigma(zz >= zD).*9.81*(rho_w_lsat-rho_w_vsat)/rho_w_lsat^2).^0.25);
alphaD_out(1) = alpha_D;
alpha0_D = linspace(0,alpha_D,length(zFDB)-length(zD_out));

alpha = zeros(1,length(zz));
alpha(1:yONB) = 0;
alpha(yONB+1:yD-1) = alpha0_D(1:end-1);
alpha(yD:end) = alphaD_out;
alphaD(yD) = alpha_D; % imposing continuity of the function

figure(7)
plot(zz,alpha,'LineWidth',1.5)
grid on
title('Void fraction profile')
xlabel('Core height [m]')
ylabel('Void fraction [-]')

%% 7 cladding inner wall temperature

% The conductivity of the cladding depends on the temperature the resolution of the conductivity integral is made analitically

aa = 1.425e-2/2;
bb = 11.45;
cc = q_vol*A_sec_f*log(D_co/D_ci)/2/pi;
beta = cc + aa*T_co_sub_JL.^2 + bb*T_co_sub_JL;

T_ci = (-bb + sqrt(bb^2 + 4*aa*beta))/2/aa;
k_cl = 11.45 + 1.425e-2*T_ci;

figure(8)
plot(T_ci, zz, LineWidth=1.5)
grid on
title('Cladding inner surface temperature profile')
xlabel('Temperature [°C]')
ylabel('Core height [m]')

%% 8-9 temperature of the fuel pellet surface and centerline
R_fact = 0.97; % Robertson factor
Ta = 20; % [°C]
P_amb_in_cl = 3e6; % [Pa] pressure of the gas at ambient temperature

% Elastic deformation of the cladding
R_ci = D_ci/2;
R_co = D_co/2;
R_p = D_p/2;
rr = linspace(0,R_p);
gap_thickness = R_ci - R_p;

T_cl_mean = (T_co_sub_JL + T_ci)/2;
    
% Cladding diameter expantion
alpha_lexp_cladd = 5.62e-6 + 3.162e-9*T_cl_mean; % [1/°C]
D_ci_exp = D_ci*alpha_lexp_cladd.*(T_cl_mean-Ta) + D_ci;
    
   
E_zr = 1.148e11 - 5.99e7*(T_cl_mean + 273.15); % [Pa] cladding Young modulus,
nu_cl = 0.43; % cladding Poisson modulus

% Ross and Stoute correlation
C_gap = 2.54e-5; % [m] coefficient taking into account fuel and cladding roughness 
AA = 0.1763e-2;
NN = 0.77163;

h_gap = 2e3; % [W/(m^2*K)]
h_c = 0;
h_rad = 0;
T_fS = T_ci + q_vol*A_sec_f./(pi*D_p*h_gap);

k_fuel_rad_ave = 4;
for ii = 1:length(zz)
   T_f_matrix(ii,:) = q_vol(ii)*R_fact./(4*k_fuel_rad_ave).*(D_p^2/4 - rr.^2) + T_fS(ii);
end
T_fC = T_f_matrix(:,1);

itermax = 1000;

iter2 = 0;
toll2 = 1e-6;
err2 = toll2*10;
while err2 > toll2 && iter2 < itermax
    iter2 = iter2+1;
    
    for ii = 1:length(zz)
            T_f_radial = T_f_matrix(ii,:);
            k_fuel = (1./(11.8 + 0.0238*T_f_radial) + 8.775e-13*T_f_radial.^3)*1e2; % [W/m*°C]
            k_fuel_rad_ave = abs(trapz(T_f_radial,k_fuel))./(T_f_radial(1)-T_fS(ii));
            T_f_matrix(ii,:) = q_vol(ii)*R_fact./(4*k_fuel_rad_ave).*(D_p^2/4 - rr.^2) + T_fS(ii);
            T_f_mean(ii) = 2/(R_p^2) * (trapz(rr,T_f_radial.*rr));
     end
     T_fC = T_f_matrix(:,1);

    iter1 = 0;
    toll1 = 1e-6;
    err1 = toll1*10;
    while err1 > toll1 && iter1 < itermax 
        iter1 = iter1+1;

        T_gap_mean = (T_fS + T_ci)/2; 
        
        % Gap thickness
        % Fuel diameter expantion
        alpha_lexp_fuel = 7.87e-6 + 3.9e-9*T_f_mean; % [1/°C]
        D_p_exp = D_p*alpha_lexp_fuel.*(T_f_mean-Ta) + D_p;

        gamma_sq = (R_co/R_ci)^2;
        P_T_in_cl = P_amb_in_cl*(T_gap_mean+273.15)/(Ta+273.15);
        delta_R_ci = R_ci./(E_zr*(gamma_sq-1)).*(P_T_in_cl*((1-nu_cl) + (1+nu_cl)*gamma_sq) - 2*gamma_sq*P_bar*1e5); 

        % Equivalent cladding diameter considering expantion+deformation
        D_ci_eq = D_ci_exp + 2*delta_R_ci;

        gap_thickness = D_ci_eq - D_p_exp;
        flag = 0;
        h_c = 0;
        negative = gap_thickness < 0;
        if any(negative)
            fprintf('fuel-cladding contact on %d nods, iteration %d\n', sum(negative), iter1);
            gap_thickness(negative) = 0;

            % Dean correlation for Contact heat transfer coefficient
            R_ci_contact = D_ci_eq(negative)/2;
            R_p_contact = D_p_exp(negative)/2;
            b_coeff = 1.36e-3;
            gamma_sq = (R_co./R_ci_contact).^2;
            U_cl = R_ci_contact./E_zr(negative).*((gamma_sq+1)./(gamma_sq-1)+nu_cl);

            % Elastic deformation of the pellet
            porosity = 0.045;
            nu_f = 0.316;
            E_f = 22.9e10 - 20.1e6*(T_f_mean+273.15) - 58.7e10*porosity;
            U_f = R_p_contact./E_f(negative)*(1+nu_f);
            p_contact = (R_p_contact-R_ci_contact)./(U_f+U_cl); %[Pa]
            
            % Contact heat transfer coefficient
            h_c = zeros(1,length(zz));
            h_c(negative) = b_coeff*p_contact;
            h_c_eff = h_c;

        end
            
        % Ross and Stoute correlation
        k_gas = AA*(T_gap_mean+273.15).^NN;
        h_gas = k_gas./(C_gap+gap_thickness);

        % Radiative heat transfer coefficient
        sigma_boltz = 5.670374419e-8;
        T_fS_K = T_fS + 273.15;
        T_ci_K = T_ci + 273.15;
        h_rad = sigma_boltz*(T_fS_K.^4 - T_ci_K.^4)./(T_fS_K - T_ci_K);

        % Total heat transfer coefficient
        h_gap = h_gas + h_c + h_rad;
    
        T_fS_new = T_ci + q_vol*A_sec_f./(pi*D_p*h_gap);
        err1 = norm(T_fS_new - T_fS)/norm(T_fS_new-T_ci);
        T_fS = T_fS_new;  

    end

    for ii = 1:length(zz)
        T_f_radial = T_f_matrix(ii,:);
        k_fuel = (1./(11.8 + 0.0238*T_f_radial) + 8.775e-13*T_f_radial.^3)*1e2; % [W/m*°C]
        k_fuel_rad_ave = abs(trapz(T_f_radial,k_fuel))./(T_f_radial(1)-T_fS(ii));
        T_f_matrix(ii,:) = q_vol(ii)*R_fact./(4*k_fuel_rad_ave).*(D_p^2/4 - rr.^2) + T_fS(ii);
    end
    
    T_fC_new = T_f_matrix(:,1);

    err2 = norm(T_fC_new - T_fC)/norm(T_fC_new'-T_ci);
    T_fC = T_fC_new;
end

T_fC_max = max(T_fC);
T_UO2_melt = 2800; %°C
fprintf('Maximum fuel centerline temperature: %.1f °C\n',T_fC_max)
if T_fC_max < T_UO2_melt
   disp("Design thermal limit respected: fuel melting point did not reach")
else
    disp("Thermal limit violated! Risk of fuel melting")
end

figure(9)
plot(T_fS,zz,LineWidth=1.5)
grid on
title('Pellet surface temperature profile')
xlabel('Temperature [°C]')
ylabel('Core height [m]')

figure(10)
plot(T_fC,zz,LineWidth=1.5)
grid on
title('Pellet centerline temperature profile')
xlabel('Temperature [°C]')
ylabel('Core height [m]')

% All axial temperture profile
figure(11)
plot(Tw,zz,LineWidth=1.5)
hold on
grid on
plot(T_co_sub_JL,zz,LineWidth=1.5)
plot(T_ci, zz, LineWidth=1.5)
plot(T_fS,zz,LineWidth=1.5)
plot(T_fC,zz,LineWidth=1.5)
ylim([-H_active/2, H_active/2])
title('Axial temperature profiles')
xlabel('Temperature [°C]')
ylabel('Core height [m]')
legend('Coolant','Outer cladding','Inner cladding','Fuel surface','Fuel centerline','Location','northeast')

%% 10 Critical heat flux

G_subch_us=G_ave*737.2; % To have in [lb/(ft^2*hr)]
G_tong=G_subch_us/1e6;
P_Mpa = P_bar*0.1;
H_active_ft=168/12;
grid_space= 168/9;
k_s=0.0704; % By interpolation 

% Uniform
F_space_grid = (P_psia/225.896)^0.5.*(1.445-0.0371.*H_active_ft).*(exp((x_eq+0.2).^2)-0.73)+k_s.*G_subch_us./1e6.*(0.038/0.019).^0.35;

q_c_15x15 = (((2.022-0.06238*P_Mpa)+(0.1722-0.01427*P_Mpa).*exp((18.177-0.5987*P_Mpa).*x_eq)).* ...
  ((0.1484-1.596.*x_eq+0.1729.*x_eq.*abs(x_eq)).*2.326.*G_ave+3271).* ...
  (1.157-0.869.*x_eq).*(0.2664+0.8357.*exp(-124.1.*D_eq)).*(0.8285+0.0003413.*(i_lsat-i_in)))*1e3;

q_c_17x17 = q_c_15x15*0.88;

q_c_15x15_factor = q_c_15x15.*F_space_grid;
q_c_17x17_factor = q_c_17x17.*F_space_grid;
qc_EU=q_c_17x17_factor(yONB:end);

% Non uniform
z_prime = linspace(0,zz(end)-zONB,length(zz(yONB:end)));
z_rel = z_prime + zONB;
FF = zeros(1,length(z_prime));
FF(1) = 1;
toll3 = 1e-10;

for tt = 2:length(z_prime)
    L_c_NU = z_prime(tt);
    y_LcNU = find(abs(zz - z_rel(tt))<toll3);
    x_LcNU = x_eq(y_LcNU);
    CC = 0.15*(1-x_LcNU)^4.31/(G_tong^0.478);
    CC_m = CC/0.0254;
    
    q_local = q_flux(y_LcNU);
    q_prime = q_flux(yONB:y_LcNU);  
    kernel = exp(-CC_m*(L_c_NU-z_prime(1:tt))); % exp[-C(L-z')]
    numerator = CC_m*trapz(z_prime(1:tt),q_prime.*kernel);
    denominator = q_local*(1-exp(-CC_m*L_c_NU));
    
    FF(tt) = numerator./denominator;
    
end

qc_NU = qc_EU./FF;

figure(12)
plot(zz,q_c_15x15_factor*1e-6,'b-','LineWidth', 2)
hold on
plot(zz,q_c_17x17_factor*1e-6,'r--','LineWidth', 2)
hold on
plot(z_rel,qc_NU*1e-6,'k-.','LineWidth',2)
hold on
plot(zz, q_flux*1e-6,'LineWidth',2)
grid on
xlabel('Core height [m]')
ylabel('Heat flux [MW/m^2]')
title('Critical Heat Flux')
legend('W3 uniform 15×15','W3 uniform 17×17','Tong non‐uniform','Heat flux','Location','Best')


%% 11 Evaluation of DNBR and minimum DNBR

DNBR = qc_NU./q_prime;

[DNBR_min,idx_min] = min(DNBR);
z_min = z_rel(idx_min);

figure(13)
plot(z_rel,DNBR,'b-','LineWidth',2)
hold on;
plot(z_min,DNBR_min,'ro','MarkerSize',8,'LineWidth',2)
yline(1.85,'--','threshold: 1.85', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right')
xlabel('Core height [m]')
ylabel('DNBR [-]')
title('Axial DNBR profile')
grid on
legend('DNBR profile','Minimum DNBR','Location','best')

if DNBR_min >= 1.85
    disp("Design thermal limit is respected (MDNBR ≥ 1.85)")
else
    fprintf('MDNBR: %.2f \n',DNBR_min)
    fprintf('Most critical point: %.2f m\n',z_min+H_active/2)
    disp("Thermal limit violated! Risk of DNB (MDNBR < 1.85)")
end
