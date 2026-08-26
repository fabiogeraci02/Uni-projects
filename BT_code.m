clc
clear
close all
%plot delle proprietà dei materiali
%conducibilità del Pyrogel XT-E
T_data = [273 373 473 573 673 773 873]; %K
kk_pyr_data = [0.02 0.023 0.028 0.035 0.046 0.064 0.089]; %W/m*K
TT_prova = linspace(273,873);
kk_pyr_interp = spline(T_data,kk_pyr_data,TT_prova);

figure()
plot(TT_prova,kk_pyr_interp,'LineWidth',2,'Color','b')
hold on
plot(T_data,kk_pyr_data,'o','LineWidth',2,'Color',[0,0,0.7])
ylim([0 0.1])
xlim([200 900])
grid on
title('Conducibilità termica Pyrogel XT-E')
xlabel('Temperatura [K]')
ylabel('k [W/m*K]')

%conducibilità dell'AISI 316
T_data = [300 400 600 800 1000]; %K
kk_steel_data = [13.4 15.2 18.3 21.3 24.2]; %W/m*K
TT_prova = linspace(300,1000);
kk_steel_interp = spline(T_data,kk_steel_data,TT_prova);

figure()
plot(TT_prova,kk_steel_interp,'LineWidth',2,'Color','r')
hold on
plot(T_data,kk_steel_data,'o','LineWidth',2,'Color',[0.8,0,0.5]) 
grid on
xlim([200 1100])
xlabel('Temperatura [K]')
ylabel('k [W/m*K]')
title('Conducibilità termica acciaio AISI 316')

%calore specifico dell'AISI 316
T_data = [300 400 600 800 1000]; %K
cp_steel_data = [468 504 550 576 602]; %W/m*K
TT_prova = linspace(300,1000);
cp_steel_interp = spline(T_data,cp_steel_data,TT_prova);

figure()
plot(TT_prova,cp_steel_interp,'LineWidth',2,'Color','r')
hold on
plot(T_data,cp_steel_data,'o','LineWidth',2,'Color',[0.8,0,0.5]) 
grid on
xlim([200 1100])
xlabel('Temperatura [K]')
ylabel('cp [J/kg*k]')
title('Calore specifico acciaio AISI 316')


pp = 11*1e5; %Pa
TT = linspace(280,910,200);
rho_He = zeros(length(TT),1);
k_He = zeros(length(TT),1);
for i = 1:length(TT)
    rho_He(i) = dhe(pp,TT(i));
    k_He(i) = condhe(pp,TT(i));
end
%Conducibilità dell'elio
figure()
plot(TT,k_He,'LineWidth',2,'Color','g')
grid on
xlabel('Temperatura [K]')
ylabel('k [W/m*K]')
title('Conducibilità termica elio')
xlim([260 1000])

%densità dell'elio
figure()
plot(TT,rho_He,'LineWidth',2,'Color',[0,0.6,0.2])
grid on
xlabel('Temperatura [K]')
ylabel('rho [kg/m^3]')
title('Densità elio')
xlim([260 1000])


%% Modello 1D-R stazionario spessore di isolante
clc
%lunghezza totale dello scambiatore + isolante
LL_is = 1.906; %m

%dicretizzazione dello spessore di isolante
r4 = 0.010; %m - raggio interno dell'isolante
r5 = 0.068; %m - raggio esterno dell'isolante
ss_is = 0.058; %m - spessore di isolante
dr = 2e-5;

rr_is = (r4:dr:r5)';
nn_is = length(rr_is);

%iterazione dovuta alla dipendenza di h_air e k_pyrogel dalla temperatura

TT_start = 903.15; %K - temperatura dello scambiatore
TT_air = 293.15; %K - temperatura dell'aria imperturbata

TT_guess = TT_start*ones(nn_is,1);

toll = 1e-8;
maxiter = 1000;
err1 = 10*toll;
iter = 0;

T_old1 = TT_guess;

while err1 > toll && iter < maxiter
    iter = iter+1;
    
    hh = myconv_fun(T_old1(end));
    kk_pyr = mycond_pyr_mod1(T_old1); 

    %creazione della matrice dei coefficienti
    subdiag1 = [-(kk_pyr(3:end)-kk_pyr(1:end-2))/4 + kk_pyr(2:end-1).*(1-dr./(2*rr_is(2:end-1)));0;0];
    supdiag1 = [0;0;(kk_pyr(3:end)-kk_pyr(1:end-2))/4 + kk_pyr(2:end-1).*(1+dr./(2*rr_is(2:end-1)))];
    main1 = [0;-2*kk_pyr(2:end-1);0];
    AA1 = spdiags([subdiag1 main1 supdiag1],-1:1,nn_is,nn_is);

    %vettore dei termini noti
    bb1 = zeros(nn_is,1);

    %condizioni al contorno
    %nodo 1 - condizione di Dirichlet (temperatura iniziale di 630 °C)
    AA1(1,1) = 1;
    bb1(1) = TT_start;

    %nodo end - condizione di Robin (scambio convettivo con aria)
    AA1(end,end) = 1 + hh*dr/kk_pyr(end);
    AA1(end,end-1) = -1;
     
    bb1(end) = hh*dr/kk_pyr(end)*TT_air;

    %risoluzione
    TT_mod1 = AA1\bb1;

    %calcolo l'errore relativo rispetto alla differenza di temperatura con
    %l'aria
    err1 = norm(TT_mod1-T_old1)/norm(TT_mod1-TT_air);
    
    %aggiorno la temperatura per il calcolo delle proprietà non costanti
    T_old1 = TT_mod1;
    
end    

figure()
plot(rr_is*1e3,TT_mod1,'LineWidth',2)
grid on
xlabel('raggio[mm]')
ylabel('Temperatura[K]')
xlim([0 70])
title('Profilo di temperatura nel Pyrogel XT-E')

figure()
plot(TT_mod1,kk_pyr,'LineWidth',2)
xlabel('Temperatura [K]')
ylabel('k pyrogel [W/m*K]')
title('Conducibilità del Pyrogel XT-E in funzione della Temperatura')
grid on

%bilancio energetico
P_in = kk_pyr(1)*(TT_mod1(1)-TT_mod1(2))/dr*pi*r4*LL_is; %W
P_out = hh*(TT_mod1(end)-TT_air)*pi*r5*LL_is; %W

%% studio di convergenza spaziale del modello 1
clc
%dicretizzazione dello spessore di isolante
r4 = 0.010; %m - raggio interno dell'isolante
r5 = 0.068; %m - raggio esterno dell'isolante

TT_start = 903.15; %K - temperatura dello scambiatore
TT_air = 293.15; %K - temperatura dell'aria imperturbata


%creo il vettore dr per le diverse discretizzazioni
vect_dr1 = sort([logspace(-5,-3,3), 2*logspace(-5,-3,3), 5*logspace(-5,-3,3)]);
err_rel = zeros(length(vect_dr1),1);

%creo un ciclo for per calcolare la soluzione variando il dr
for jj = 1:length(vect_dr1)
    
    dr = vect_dr1(jj);
    rr_is = (r4:dr:r5)';
    nn_is = length(rr_is);
    
    TT_guess = TT_start*ones(nn_is,1);

    toll = 1e-8;
    maxiter = 1000;
    err1 = 10*toll;
    iter = 0;
    
    T_old1 = TT_guess;

    while err1 > toll && iter < maxiter
        iter = iter+1;
    
        hh = myconv_fun(T_old1(end));
        kk_pyr = mycond_pyr_mod1(T_old1); 
    
    
        %creazione della matrice dei coefficienti
        subdiag1 = [-(kk_pyr(3:end)-kk_pyr(1:end-2))/4 + kk_pyr(2:end-1).*(1-dr./(2*rr_is(2:end-1)));0;0];
        supdiag1 = [0;0;(kk_pyr(3:end)-kk_pyr(1:end-2))/4 + kk_pyr(2:end-1).*(1+dr./(2*rr_is(2:end-1)))];
        main1 = [0;-2*kk_pyr(2:end-1);0];
        AA1 = spdiags([subdiag1 main1 supdiag1],-1:1,nn_is,nn_is);
    
        %vettore dei termini noti
        bb1 = zeros(nn_is,1);
    
        %condizioni al contorno
        %nodo 1 - condizione di Dirichlet (temperatura iniziale di 630 °C)
        AA1(1,1) = 1;
        bb1(1) = TT_start;
    
        %nodo end - condizione di Robin (scambio convettivo con aria)
        AA1(end,end) = 1 + hh*dr/kk_pyr(end);
        AA1(end,end-1) = -1;
         
        bb1(end) = hh*dr/kk_pyr(end)*TT_air;
    
        %risoluzione
        TT_mod1_conv = AA1\bb1;
          %calcolo l'errore relativo rispetto alla differenza di temperatura con
        %l'aria
        err1 = norm(TT_mod1_conv-T_old1)/norm(TT_mod1_conv-TT_air);
        
        %aggiorno la temperatura per il calcolo delle proprietà non costanti
        T_old1 = TT_mod1_conv;
        
    end    

    if jj == 1
        dr_rif = dr;
        TT_rif = TT_mod1_conv;
    end
    %compongo il vettore degli errori di discretizzazione
    err_rel(jj) = norm(TT_mod1_conv-TT_rif(1:round(dr/dr_rif):end))/norm(TT_rif(1:round(dr/dr_rif):end)-TT_air);

    K_AA = condest(AA1);
    err_arr(jj) = 2*eps*K_AA/(1-K_AA*eps);

end
figure()
loglog(vect_dr1,err_rel,'o-r')
hold on
loglog(vect_dr1,err_arr,'d--b')
grid on
xlabel('\Deltar [m]')
ylabel('errore relativo [-]')
title('Convergenza spaziale')
legend('err. di troncamento','Up.Bound err.arrotondamento','Location','northwest')


%% Modello 1D-R transitorio sul sistema isolante + baionetta
clc
%discretizzazione del sistema
dr1 = 2e-6; 
dr2 = 2e-5;

r1 = 5.5e-3; %m
r2 = 7e-3;
r3 = 7.5e-3;
r4 = 10e-3;
r5 = 68e-3;

vect1 = (r1:dr1:r2)';
nn1 = length(vect1);

vect2 = (r2+dr1:dr1:r3)';
nn2 = length(vect2)+nn1;

vect3 = (r3+dr1:dr1:r4)';
nn3 = length(vect3)+nn2;

vect4 = (r4+dr2:dr2:r5)';
nn4 = length(vect4)+nn3;

rr = [vect1;vect2;vect3;vect4];
nn = length(rr);

%Creo la condizione iniziale del transitorio con il vettore temperatura T
%uniforme nella baionetta + T isolante

TT_start = 903.15; %K - temperatura dello scambiatore
TT_air = 293.15; %K - temperatura dell'aria imperturbata

T0_is = TT_mod1; %salvo la soluzione del primo modello
T0_baionetta = TT_start*ones(nn3-1,1);

%Temperatura iniziale
T0 = [T0_baionetta; T0_is];

%definisco tutte le proprietà dei materiali INDIPENDENTI dalla temperatura
%acciaio AISI 316/316L
rho_steel = 8238; %kg/m^3

%elio stagnante (calcolate a 1 atm)
cp_He = 5193; %J/kg*K

%isolante Pyrogel XT-E
rho_pyr = 200.23; %kg/m^3
cp_pyr = 2300; %J/kg*K

%discretizzazione temporale
dt = 60;

%vettore dei tempi in cui vengono stampati i profili di temperatura
vect_time1 = [dt 60*dt 180*dt 360*dt 600*dt 1200*dt 124860];
i = 1;

%ciclo fin quando non arrivo a un profilo di temperatura asintotico
%(stazionario finale)
toll = 1e-5;
err2 = 10*toll;
iter = 0;
time = 0;
maxiter = 15000;
Told2 = T0;

while err2>toll && iter<maxiter
   
    iter = iter+1;
    time = time+dt;

    %definisco tutte le proprietà dei materiali DIPENDENTI dalla
    %temperatura sfruttando il metodo del "Frozen coefficient"
    %acciaio AISI 316/316L
    [k_steel1,k_steel2] = mycond_steel(Told2);
    [cp_steel1,cp_steel2] = mycapacity_steel(Told2);

    %elio stagnante (calcolate a 1 atm)
    [k_He,rho_He] = myproperties_he(Told2(nn1+1:nn2));

    %isolante Pyrogel XT-E
    k_pyr = mycond_pyr_mod2(Told2);

    %ricavo il coefficiente di scambio convettivo
    hh = myconv_fun(Told2(end));

    %Creo il vettore della diffusività termica di tutto il componente
    %usando i vettori logici
    %conducibilità termica
    kk = k_steel1.*(rr<=r2) + k_steel2.*(rr>r3 & rr<=r4) + k_pyr.*(rr>r4);
    kk = [kk(1:nn1);k_He;kk(nn2+1:end)];
    
    %calore specifico
    cp = cp_steel1.*(rr<=r2) + cp_He*(rr>r2 & rr<=r3) + cp_steel2.*(rr>r3 & rr<=r4) + cp_pyr*(rr>r4);
    
    %densità massica
    rho = rho_steel*(rr<=r2) + rho_steel*(rr>r3 & rr<=r4) + rho_pyr*(rr>r4);
    rho = [rho(1:nn1);rho_He;rho(nn2+1:end)];
    
    %diffusività termica
    alpha = kk./(cp.*rho); 

    dr_logic = dr1*(rr<=r4) + dr2*(rr>r4);
    aa = (alpha./(dr_logic.^2))*dt;
    
    %Definisco la matrice dei coefficienti e il vettore dei termini noti
    subdiag2 = [-aa(2:end-1).*((1-dr_logic(2:end-1)./(2*rr(2:end-1)))); 0; 0];
    supdiag2 = [0; 0; -aa(2:end-1).*((1+dr_logic(2:end-1)./(2*rr(2:end-1))))];
    main2 = [0; 1+2*aa(2:end-1); 0];

    AA2 = spdiags([subdiag2 main2 supdiag2],-1:1,nn,nn);
    bb2 = Told2;

    %condizioni al contorno
    %nodo 1 - Condizione di adiabaticità(Neumann)
    AA2(1,1) = 1;
    AA2(1,2) = -1;
    bb2(1) = 0;

    %nodo end - Condizione di scambio termico convettivo(Robin)
    AA2(end,end) = 1 + hh*dr_logic(end)/k_pyr(end);
    AA2(end,end-1) = -1;
    bb2(end) = (hh*dr_logic(end)/k_pyr(end))*TT_air;

    %condizioni all'interfaccia (Continuità del flusso)
    %interfaccia 1-2 (acciaio-elio)
    AA2(nn1,nn1) = 1 + k_steel1(end)/k_He(1);
    AA2(nn1,nn1-1) = - k_steel1(end)/k_He(1);
    AA2(nn1,nn1+1) = -1;
    bb2(nn1) = 0;

    %interfaccia 2-3 (elio-acciaio)
    AA2(nn2,nn2) = 1 + k_He(end)/k_steel2(1);
    AA2(nn2,nn2-1) = -  k_He(end)/k_steel2(1);
    AA2(nn2,nn2+1) = -1;
    bb2(nn2) = 0;

    %interfaccia 3-4 (acciaio-isolante)
    AA2(nn3,nn3) = 1 + (k_steel2(end)/k_pyr(1))*(dr2/dr1);
    AA2(nn3,nn3-1) = - (k_steel2(end)/k_pyr(1))*(dr2/dr1);
    AA2(nn3,nn3+1) = -1;
    bb2(nn3) = 0;

    %soluzione
    TT_mod2 = AA2\bb2;

    %errore relativo
    err2 = norm(TT_mod2-Told2)/norm(TT_mod2-TT_start);
    
    %aggiorno il vettore temperatura al dt precedente
    Told2 = TT_mod2;
    
    %stampo solo alcuni profili di temperatura durante il transitorio
    if time == vect_time1(i)
        figure(9)
        plot(rr*1e3,Told2,LineWidth=1.5,DisplayName=strcat('time= ',num2str(round(vect_time1(i)/3600,2)),'h'));
        hold on
        grid on
        xlabel('raggio [mm]')
        ylabel('Temperatura [K]')
        title('Profili di temperatura durante il transitorio')
        legend('-dynamiclegend')
        i = i+1;
    end
    
    %Calcolo della potenza
    Potenza_dissipata(iter) = hh*(Told2(end)-TT_air)*pi*r5*LL_is; %W

    %creo un vettore che memorizzi il  valore del coefficiente di scambio
    %convettivo durante il transitorio
    vect_hh(iter) = hh; 

end

%grafico della potenza dissipata
vect_time2 = dt:dt:time;
figure()
plot(vect_time2/3600,Potenza_dissipata,LineWidth=2)
xlabel('time [h]')
ylabel('Potenza dissipata [W]')
title('Potenza dissipata durante il transitorio')
grid on

%visualizzazione del profilo nella sola baionetta
figure()
plot([vect1;vect2;vect3]*1e3,TT_mod2(1:nn3),LineWidth=2)
grid on
xlabel('raggio [mm]')
ylabel('Temperatura [K]')
title('Profilo di temperatura finale della baionetta')

%grafico della diffusività termica nella baionetta
figure()
semilogy([vect1;vect2;vect3]*1e3,alpha(1:nn3),'LineWidth',2)
ylim([1e-7 1e-3])
grid on

%andamento del valore di h
figure()
plot(vect_time2/3600,vect_hh,'LineWidth',2)
grid on
xlabel('time [h]')
ylabel('h [W/m^2*K]')
title('Evoluzione del coeffciente h')


%% Definizione delle funzioni per le proprietà dei materiali dipendenti dalla Temperatura
function[k_pyr] = mycond_pyr_mod1(TT)
    T_data1 = [273 373 473 573 673 773 873]; %K
    kk_pyr_data = [0.02 0.023 0.028 0.035 0.046 0.064 0.089]; %W/m*K
    k_pyr = spline(T_data1,kk_pyr_data,TT);
end

function[k_pyr] = mycond_pyr_mod2(TT)
    T_data1 = [273 373 473 573 673 773 873]; %K
    kk_pyr_data = [0.02 0.023 0.028 0.035 0.046 0.064 0.089]; %W/m*K
    k_pyr = spline(T_data1,kk_pyr_data,TT);
end

function[k_steel1,k_steel2] = mycond_steel(TT)
    T_data2 = [300 400 600 800 1000]; %K
    kk_steel_data = [13.4 15.2 18.3 21.3 24.2]; %W/m*K
    k_steel1 = spline(T_data2,kk_steel_data,TT);
    k_steel2 = spline(T_data2,kk_steel_data,TT);
end

function[cp_steel1,cp_steel2] = mycapacity_steel(TT)
    T_data2 = [300 400 600 800 1000]; %K
    cp_steel_data = [468 504 550 576 602]; %J/kg*K
    cp_steel1 = spline(T_data2,cp_steel_data,TT);
    cp_steel2 = spline(T_data2,cp_steel_data,TT);
end

function[k_He,rho_He] = myproperties_he(TT)
    pp = 11*1e5; %Pa
    rho_He = zeros(length(TT),1);
    k_He = zeros(length(TT),1);
    for i = 1:length(TT)
        rho_He(i) = dhe(pp,TT(i));
        k_He(i) = condhe(pp,TT(i));
        %cp_He(i) = cphe(pp,TT(i));
    end
end

function[hh] = myconv_fun(Tguess)
    %proprietà dell'aria calcolate a 20°C e 1 atm
    T_out = 293.15; %K
    kk_air = 0.0259; %W/m*K
    mu = 18.2e-6; %Pa*s
    cp_air = 1006; %J/kg*K
    rho_air = 101325/(287*T_out); %kg/m^3
    alpha_air = kk_air/(rho_air * cp_air); %m^2/s
    ni = mu/rho_air; %m^2/s
    beta = 1/T_out; %1/K

    %Dati geometrici
    LL = 1.906; %m

    %valuto il numero di Rayleigh per scegliere la correlazione corretta
    Ra_L = 9.81*beta*(Tguess-T_out)*LL^3/(ni*alpha_air);

    if Ra_L >= 1e7 && Ra_L <= 1e9
        Nu_L = 0.67*Ra_L^0.25;
    elseif Ra_L > 1e9
        Nu_L = 0.0782*Ra_L^0.357;
    end

    %trovato il numero di Nusselt ricavo il coefficiente di scambio
    %convettivo

    hh = Nu_L*kk_air/LL; %W/m^2*K 
end