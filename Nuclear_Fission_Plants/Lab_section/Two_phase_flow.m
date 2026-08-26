clear all;
close all;
clc;
set(0,'defaultaxesfontsize',12);
% set(0,'defaultlinelinewidth',3);
set(0,'defaultlegendlocation','best');

%% 1) Upload data
Datas = readtable('Phase flow datas.xlsx','VariableNamingRule','preserve');
Test_num               = Datas{:,1};
Water_diaphragm        = Datas{:,2};
Air_flow_rate          = Datas{:,3};
Rotameter_air_pressure = Datas{:,4};
Water_pressure_drop    = Datas{:,5};
Relative_pressure      = Datas{:,6};
Diff_pressure_type     = Datas{:,7};
Pipe_pressure_drop     = Datas{:,8};
Water_temperature      = Datas{:,9};
Water_mass             = Datas{:,10};
Flow_pattern           = Datas{:,11};

%% 2) Parameter tube e fluid
D_in  = 0.026;
Area  = pi*(D_in/2)^2;
rho_L = 1000;
rho_G = ((Relative_pressure+1.013)*1e5*0.029)./(8.314*(Water_temperature+273.15));
sigma = 0.072;
g     = 9.81;

%% 3) gas flow and superficial velocity
rho_ref = 1.204;                             % kg/m³ a 20°C e 1 atm
Qs      = str2double(Air_flow_rate);        % Nm³/h
pe_rel  = str2double(Rotameter_air_pressure);% bar
W_e     = (rho_ref .* sqrt((pe_rel+1.013)/1.013) .* Qs)/3600;  % kg/s
V_air   = W_e ./ rho_G;                     % m³/s
jG      = V_air ./ Area;                    % m/s

%% 4) Calculation of liquid flow rate and customized beta
N         = length(Water_diaphragm);
m_L       = nan(N,1);
jL_full   = nan(N,1);
beta_vec  = nan(N,1);
rho_w     = 1000;
jL_full = ones(14,1);
for i = 1:N
    
    dpPa = Water_pressure_drop(i) * 100;  % mbar → Pa
    if strcmpi(Water_diaphragm{i}, 'S')
        d = 0.0084; a = 0.650729; B = 0.01415274; beta_val = 0.5;
    elseif strcmpi(Water_diaphragm{i}, 'M')
        d = 0.0153; a = 0.728193; B = 0.00375718; beta_val = 0.7217;
    else
        error('Tipo diaframma non valido in riga %d', i);
    end
    A_dia        = pi*(d)^2/4;
    m_L(i)       = a * A_dia * sqrt(2 * rho_w * dpPa) + B;
    jL_full(i)   = m_L(i) / (rho_L * Area);
    beta_vec(i)  = beta_val;
end

%% 5) Void fraction sperimental
M_res     = 0.116;
M_l0      = 1.014;
M_l       = Water_mass./1000+M_res;
alpha_exp = 1-(M_l./M_l0);
m_G       = rho_G.*V_air;
xx        = m_G./(m_G+m_L);

%% 6) 
maskS    = strcmpi(Water_diaphragm,'S');
maskM    = strcmpi(Water_diaphragm,'M');
mask     = maskS | maskM;
isSlug   = contains(Flow_pattern,'Slug','IgnoreCase',true);
isChurn  = contains(Flow_pattern,'Churn','IgnoreCase',true);
isAnnular= contains(Flow_pattern,'Annular','IgnoreCase',true);

%% 7) Coordinates
num = (jG .* (rho_L)^(1/2));
den = (g*(rho_L-rho_G)*sigma).^(0.25);
x_TD = num./den;
y_TD      = jL_full./jG;
beta_SC = jL_full./(jL_full+jG);   
x_TD_sc   = jG./sqrt(g*D_in);
y_TD_sc   = beta_SC;
mu_L      = 1e-3; mu_G = 1.8e-5;
Re_L      = rho_L.*jL_full.*D_in./mu_L;
Re_G      = rho_G.*jG.*D_in./mu_G;
f_L       = 0.079./(Re_L.^0.25);
f_G       = 0.079./(Re_G.^0.25);
dpdz_L    = rho_L.*jL_full.^2;
dpdz_G    = rho_G.*jG.^2;
x_TD_sc_a = sqrt(dpdz_L./dpdz_G);
y_TD_sc_a = jG .*(rho_G).^(0.5) ./ (g*(rho_L-rho_G)*sigma).^(0.25);

%% HEWITT–ROBERTS MAP (log–log)

chi = rho_L .* jL_full.^2;     % (kg·s⁻²·m⁻¹)
Y   = rho_G .* jG.^2;          % (kg·s⁻²·m⁻¹)

labels = cell(size(Flow_pattern));

for k = 1:numel(Flow_pattern)
    
    fp = lower(strtrim(Flow_pattern{k}));
    fp = regexprep(fp,'\([^)]*\)','');
    fp = regexprep(fp,'[\/\\,&]','-');
    fp = regexprep(fp,'\s*-+\s*','-');
    fp = strrep(fp,' ','');

    initials = {};

    if contains(fp,'bubble') || contains(fp,'bubbly')
        initials{end+1} = 'B';
    end
    if contains(fp,'slug')
        initials{end+1} = 'S';
    end
    if contains(fp,'churn')
        initials{end+1} = 'C';
    end
    if contains(fp,'annular')
        initials{end+1} = 'A';
    end
    if ~isempty(initials)
        labels{k} = strjoin(initials,'-');
    else
        labels{k} = '?';         
    end
end

% plot
figure('Name','Hewitt–Roberts Map','Color','w');
hold on

loglog(chi, Y, 'o', 'MarkerEdgeColor','none','MarkerFaceColor',[.8 .8 .8], 'MarkerSize',6);

% Letters
for k = 1:numel(chi)
    text(chi(k), Y(k), labels{k},'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontSize',10);
end

set(gca, 'XScale','log', 'YScale','log','XMinorGrid','on', 'YMinorGrid','on','GridLineStyle','-', 'MinorGridLineStyle',':', ...
         'GridAlpha',0.30, 'MinorGridAlpha',0.18, 'TickDir','out', 'Box','on');

xlim([1e-1 1e6]); ylim([1e-2 1e6]);
xlabel('\chi = \rho_L j_L^2  (kg\cdots^{-2}\cdotm^{-1})');
ylabel('Y = \rho_G j_G^2  (kg\cdots^{-2}\cdotm^{-1})');
title('Hewitt–Roberts Flow-Pattern Map');

%% TAITEL–DUCKLER MAP
figure('Name','Taitel–Duckler Map','Color','w'); hold on;
semilogx(x_TD, y_TD, 'o', 'MarkerEdgeColor','none','MarkerFaceColor',[.8 .8 .8], 'MarkerSize',6);
for k = 1:numel(x_TD)
    text(x_TD(k), y_TD(k), labels{k},'HorizontalAlignment','center', ...
         'VerticalAlignment','middle','FontWeight','bold','FontSize',10);
end
set(gca,'XScale','log','YScale','linear','XMinorGrid','on','YMinorGrid','on', ...
         'GridLineStyle','-','MinorGridLineStyle',':','GridAlpha',0.3,'MinorGridAlpha',0.18, ...
         'TickDir','out','Box','on');
xlim([0.1 100]); ylim([0.1 2.5]);
xlabel('j_G\surd\rho_L/(g\Delta\rho\sigma)^{1/4}'); ylabel('j_L/j_G');
title('Taitel–Duckler map – transition between Bubbly and Slug/Churn flow');

%% SLUG→CHURN TRANSITION
figure('Name','Slug\rightarrowChurn Transition','Color','w'); hold on;
semilogx(x_TD_sc, y_TD_sc, 'o', 'MarkerEdgeColor','none','MarkerFaceColor',[.8 .8 .8], 'MarkerSize',6);
for k = 1:numel(x_TD_sc)
    text(x_TD_sc(k), y_TD_sc(k), labels{k},'HorizontalAlignment','center','VerticalAlignment','middle', ...
         'FontWeight','bold','FontSize',10);
end
set(gca,'XScale','log','YScale','linear','XMinorGrid','on','YMinorGrid','on', ...
         'GridLineStyle','-','MinorGridLineStyle',':','GridAlpha',0.3,'MinorGridAlpha',0.18,'TickDir','out','Box','on');
xlim([1 1e4]); ylim([0 1]);
xlabel('j_G/\surd(gD)'); ylabel('\beta');
title('Transition between Slug and Churn flow');

%% SLUG/CHURN→ANNULAR TRANSITION
figure('Name','Slug/Churn\rightarrowAnnular','Color','w'); hold on;
loglog(x_TD_sc_a, y_TD_sc_a, 'o', 'MarkerEdgeColor','none', ...
                     'MarkerFaceColor',[.8 .8 .8], 'MarkerSize',6);
for k = 1:numel(x_TD_sc_a)
    text(x_TD_sc_a(k), y_TD_sc_a(k), labels{k},'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontSize',10);
end
set(gca,'XScale','log','YScale','log','XMinorGrid','on','YMinorGrid','on','GridLineStyle','-','MinorGridLineStyle',':', ...
         'GridAlpha',0.3,'MinorGridAlpha',0.18,'TickDir','out','Box','on');
xlim([1e-2 1e4]); ylim([1e-3 1e3]);
xlabel('\surd(dp/dz_L/dp/dz_G)');
ylabel('j_G\surd\rho_G/(g\Delta\rho\sigma)^{1/4}');
title('Taitel–Duckler map – transition between Annular and Slug/Churn flow');


%% (2) Void fraction
alpha_hom    = 1 ./ (1+((1-xx)./xx).*(rho_G./rho_L));
S_Zivi       = (rho_L./rho_G).^(1/3);
alpha_Zivi   = 1 ./ (1+((1-xx)./xx).*(rho_G./rho_L).*S_Zivi);
S_ch         = sqrt(1-xx.*(1-rho_L./rho_G));
alpha_ch     = 1 ./ (1+((1-xx)./xx).*(rho_G./rho_L).*S_ch);
beta         = (rho_L.*xx)./(rho_L.*xx+rho_G.*(1-xx));
yy           = (1-beta)./(1+beta);
Re           = (rho_L.*jL_full.*D_in)./1e-3;
We           = (rho_L.*jL_full.^2.*D_in)./sigma;
E1           = 1.578.*Re.^(-0.19).*(rho_L./rho_G).^0.22;
E2           = 0.0273.*We.^(-0.51).*(rho_L./rho_G).^(-0.08);
S_ci         = 1+E1.*((yy./(1+yy.*E2))-yy.*E2).^0.5;
alpha_cise   = 1 ./ (1+((1-xx)./xx).*(rho_G./rho_L).*S_ci);
j_tot        = jL_full + jG;
ugj_b        = 1.41.*((sigma.*(rho_L-rho_G))./rho_L.^2).^0.25;
alpha_drift_b= jG ./ (1.13.*j_tot + ugj_b);
ugj_p        = 0.35.*sqrt(((rho_L-rho_G).*g.*D_in)./rho_L);
alpha_drift_p= jG./(1.2.*j_tot+ugj_p);

figure(5); clf; hold on; grid on; axis square;
t = linspace(0,1,100);
plot(t,t,'r-','LineWidth',1.5,'DisplayName','y = x');
plot(t,1.2*t,'r--','LineWidth',1.2,'DisplayName','+20%');
plot(t,0.8*t,'r--','LineWidth',1.2,'DisplayName','-20%');
scatter(alpha_exp,alpha_hom,75,'s','filled','DisplayName','Homogeneous');
scatter(alpha_exp,alpha_Zivi,75,'o','filled','DisplayName','Zivi');
scatter(alpha_exp,alpha_ch,75,'^','filled','DisplayName','Chisholm');
scatter(alpha_exp,alpha_cise,75,'d','filled','DisplayName','CISE');
scatter(alpha_exp,alpha_drift_b,75,'p','filled','DisplayName','Drift Bubbly');
scatter(alpha_exp,alpha_drift_p,75,'h','filled','DisplayName','Drift Plug');
xlabel('Experimental Void Fraction');
ylabel('Calculated Void Fraction');
title('Void fraction');
legend('Location','best');

%% (3) Evaluation of the experimental pressure difference along the vertical pipe

height = 1.5; %[m]
dp_AB = NaN(N,1);
msg = "Pressure drop in the test section oscillates between positive and negative";

for jj = 1:N
    type = strtrim(upper(Diff_pressure_type{jj}));
    
    % Avoid error if there is no pressure
    if isnan(Pipe_pressure_drop(jj))
        continue;
    end

    dp_transducer = Pipe_pressure_drop(jj)*100;

    if contains(type, 'D-C')
        dp_AB(jj) = (rho_L*g*height)-dp_transducer;
    elseif contains(type, 'C-D')
        dp_AB(jj) = dp_transducer+(rho_L*g*height); 
    end
end

% In order to visualize the message
dp_AB_display = strings(N,1);
for i = 1:N
    if isnan(dp_AB(i))
        dp_AB_display(i) = msg;
    else
        dp_AB_display(i) = sprintf('%.2f', dp_AB(i));
    end
end

% Table
T = table(string(Diff_pressure_type), dp_AB_display,'VariableNames', {'Pressure Type', 'DeltaP_AB'});

f = uifigure('Name','Results Δp_{AB}');
uitable(f,'Data', T,'Position', [20 20 500 400]); 

%% (4) Calculate elevation pressure changes using the void fraction evaluated 
% by the homogeneous model, Chisholm, CISE and drift flux correlations. 

% Rho average
rho_ave_hom     = alpha_hom.*rho_G+(1-alpha_hom)*rho_L;
rho_ave_ch      = alpha_ch.*rho_G +(1-alpha_ch)*rho_L;
rho_ave_cise    = alpha_cise.* rho_G+(1-alpha_cise)*rho_L;
rho_ave_drift_b = alpha_drift_b.*rho_G+(1-alpha_drift_b)*rho_L;
rho_ave_drift_p = alpha_drift_p.*rho_G+(1-alpha_drift_p)*rho_L;
rho_ave_exp = alpha_exp.*rho_G+(1-alpha_exp)*rho_L ;

% Elevation pressure changes
dp_elev_hom     = rho_ave_hom*g*height;
dp_elev_ch      = rho_ave_ch*g*height;
dp_elev_cise    = rho_ave_cise*g*height;
dp_elev_drift_b = rho_ave_drift_b*g*height;
dp_elev_drift_p = rho_ave_drift_p*g*height;
dp_elev_exp     = rho_ave_exp*g*height;

% Table
Tab_elev_pressure = table(dp_elev_exp,dp_elev_hom,dp_elev_ch,dp_elev_cise,dp_elev_drift_b,dp_elev_drift_p,'VariableNames', {'Experimental_Elevation','Elevation_Homogeneous','Elevation_Chisholm','Elevation_CISE', ...
        'Elevation_Drift_Bubbly','Elevation_Drift_Plug'});

f1 = uifigure('Name','Elevation Pressure Changes');
uitable(f1,'Data', Tab_elev_pressure,'Position', [20 20 400 400]); 


%% (5) Calculate friction pressure drops  by the homogeneous and Friedel models
dp_friction_exp = dp_AB-(rho_ave_exp*g*height);
dp_friction_hom = dp_AB-(rho_ave_hom*g*height);

Re_L = rho_L.*jL_full.*D_in./mu_L;
f_L  = 0.079./Re_L.^0.25;   
f_L(Re_L<=3000)=64./Re_L(Re_L<=3000);
dpdz_L = f_L.*rho_L.*jL_full.^2./(2*D_in);

% 2) factor φ²
Re_G = rho_G.*jG.*D_in./mu_G;
f_G  = 0.079./Re_G.^0.25;   
f_G(Re_G<=3000)=64./Re_G(Re_G<=3000);
rho_h = 1./(xx./rho_G+(1-xx)./rho_L);
E = (1-xx).^2+xx.^2.*(rho_L./rho_G).*(f_G./f_L);
F = xx.^0.78.*(1-xx).^0.224;
H = (rho_L./rho_G).^0.91.*(mu_G/mu_L).^0.19.*(1-mu_G/mu_L).^0.7;
Fr= j_tot.^2./(g*D_in*rho_h);
We= j_tot.^2.*D_in./(sigma*rho_h);
phi2 = E+3.24.*F.*H./(Fr.^0.045.*We.^0.035);

% 3)  Δp_friction Friedel
dpdz_friedel        = dpdz_L .* phi2;
dp_friction_friedel = dpdz_friedel * height;

% Table
Tab_friction_pressure = table(dp_friction_hom, dp_friction_friedel,'VariableNames', {'Friction ΔP Homogeneous','Friction ΔP Friedel'});
f2 = uifigure('Name','Friction Pressure Drops');
uitable(f2,'Data',Tab_friction_pressure,'Position',[20 20 450 400]);

%% (6) Calculate the total pressure change 
dp_total_exp = dp_elev_exp+dp_friction_exp;
% Homogeneous for elevation and homogeneous for friction
dp_total_hom = dp_elev_hom+dp_friction_hom;

% Cise for elevation and Friedel for friction
dp_total_friedel_cise = dp_elev_cise+dp_friction_friedel;

% Cise for elevation and Chisholm for friction
rho_ave_chisholm = alpha_ch.*rho_G+(1-alpha_ch).*rho_L;
dp_friction_chisholm = dp_AB-(rho_ave_chisholm*g*height);

dp_total_chisholm_cise = dp_elev_cise+dp_friction_chisholm;

% Drift Flux for elevation and Friedel for friction
dp_total_friedel_drift = dp_elev_drift_b+dp_friction_friedel;

lens = [length(dp_total_hom), length(dp_total_friedel_cise), length(dp_total_chisholm_cise), length(dp_total_friedel_drift)];
min_len = min(lens);


dp_total_hom = dp_total_hom(1:min_len);
dp_total_friedel_cise = dp_total_friedel_cise(1:min_len);
dp_total_chisholm_cise = dp_total_chisholm_cise(1:min_len);
dp_total_friedel_drift = dp_total_friedel_drift(1:min_len);


dp_total_hom = dp_total_hom(:);
dp_total_friedel_cise = dp_total_friedel_cise(:);
dp_total_chisholm_cise = dp_total_chisholm_cise(:);
dp_total_friedel_drift = dp_total_friedel_drift(:);

Tab_total_pressure = table(dp_total_hom, dp_total_friedel_cise, dp_total_chisholm_cise, dp_total_friedel_drift, ...
    'VariableNames', {'Total_Homogeneous', 'Total_Friedel_CISE', 'Total_Chisholm_CISE', 'Total_Friedel_DriftFlux'});

f3 = uifigure('Name','Total Pressure Change');
uitable(f3, 'Data', Tab_total_pressure, 'Position', [20 20 550 400]);


%% (7) Comparison: Experimental vs Calculated Total Pressure Drop – Clean Excel-like style
idx_hom           = ~isnan(dp_total_exp) & ~isnan(dp_total_hom);
idx_friedel_cise  = ~isnan(dp_total_exp) & ~isnan(dp_total_friedel_cise);
idx_chisholm_cise = ~isnan(dp_total_exp) & ~isnan(dp_total_chisholm_cise);
idx_friedel_drift = ~isnan(dp_total_exp) & ~isnan(dp_total_friedel_drift);

% max per axis
maxX = max(dp_total_exp(~isnan(dp_total_exp))) * 1.1;
all_y = [dp_total_hom, dp_total_friedel_cise, dp_total_chisholm_cise, dp_total_friedel_drift];

maxY = max(all_y(~isnan(all_y))) * 1.1;

figure; 
hold on; 
grid on; 
box on;

% Red lines
x_ref = linspace(0, maxX, 200);
plot(x_ref, x_ref,     'r-',  'LineWidth', 1.8, 'DisplayName', 'diagonal');
plot(x_ref, 1.2*x_ref, 'r--', 'LineWidth', 1.2, 'DisplayName', '+20%');
plot(x_ref, 0.8*x_ref, 'r--', 'LineWidth', 1.2, 'DisplayName', '-20%');

% Points
scatter(dp_total_exp(idx_hom),dp_total_hom(idx_hom),60, 'o', 'filled','MarkerFaceColor', [0.30 0.65 0.90], 'DisplayName', 'Homogeneous');        
scatter(dp_total_exp(idx_friedel_cise),  dp_total_friedel_cise(idx_friedel_cise), 60, 'o', 'filled','MarkerFaceColor', [0.20 0.60 0.25], 'DisplayName', 'Friedel + CISE');       
scatter(dp_total_exp(idx_chisholm_cise), dp_total_chisholm_cise(idx_chisholm_cise), 60, 'o', 'filled','MarkerFaceColor', [0.85 0.33 0.10], 'DisplayName', 'Chisholm + CISE');      
scatter(dp_total_exp(idx_friedel_drift), dp_total_friedel_drift(idx_friedel_drift), 60, 'o', 'filled','MarkerFaceColor', [0.49 0.18 0.56], 'DisplayName', 'Friedel + Drift Flux'); 

% lable and title
xlabel('Experimental pressure drop (Pa)');
ylabel('Calculated pressure drop (Pa)');
title('Comparison: Experimental vs Calculated Total Pressure Drop');

% Limits & legend
xlim([0 maxX]);
ylim([0 maxY]);
legend('Location', 'eastoutside');

