close all
clear all
clc


%% ============ Load Data ============
load 'input_data/covidbahia_test_Vac'
load 'input_data/Psi_Ba13set'

npd = length(Psi);
Nfo = 8;            % Filter order
for k = 1:npd
    if (k-Nfo) < 0
        Psif(k) = 1/k*sum(Psi(1:k));
    else
        Psif(k) = 1/Nfo*sum(Psi(k-Nfo+1:k));
    end
end

coviddata=covidbahia;

Cases_dado      = coviddata(:,1)';
R_dado          = coviddata(:,7)';
D_dado          = coviddata(:,5)';
H_dado          = coviddata(:,2)';
U_dado          = coviddata(:,3)';
V_dado          = coviddata(:,8)';
data_in = datetime(2020,03,06);
Nw_cases        = diff(Cases_dado);
Nw_deaths       = diff(D_dado);

t_vac           = 100; %% day started the vaccination
%% ============ Initial Condition ============

D0          = 0;
N0          = 14930634;
R0          = 0;
H0          = 0;
HV0         = 0;
U0          = 0;
UV0         = 0;
V0          = 1000/N0;  %% new compartment
Rv0         = 0; %% new compartment
Is0         = 2.015439771376298e-06;
Ia0         = 1.8028646508967777e-06;
Iav0        = 1.8028646508967777e-06; %% new compartment
Isv0        = 1.8028646508967777e-06; %% new compartment
E0          = 1.7639153732952095e-06;
Ev0         = 1.7639153732952095e-06; %% new compartment
S0          = (1-Is0-Ia0-E0);
Nw0         = 0;
NwV0        = 0;
Model_0     = [S0,E0,V0,Ev0,Ia0,Is0,Iav0,Isv0,H0,HV0,U0,UV0,R0,Rv0,D0,Nw0,NwV0];

%% ============ Model ============
    
tsim = length(U_dado)-1; 
    
[t,y] = ode45(@(t,y)covid_odes(t,y,Psif),0:1:tsim,Model_0');

Sm=N0.*y(:,1);
Em=N0.*y(:,2);
Vm=N0.*y(:,3);
Evm=N0.*y(:,4);
Iam=N0.*y(:,5);
Ism=N0.*y(:,6);
Iavm=N0.*y(:,7);
Isvm=N0.*y(:,8);
Hm=N0.*y(:,9);
Hv=N0.*y(:,10);
Um=N0.*y(:,11);
Uv=N0.*y(:,12);
Rm=N0.*y(:,13);
Rvm=N0.*y(:,14);
Dm=N0.*y(:,15);
Nwm=N0.*y(:,16);
NwVm=N0.*y(:,17);

%% ============ DATA ============

t_p = data_in + caldays(0:(length(U_dado)-1));
t_dado=1:1:length(H_dado);
t_dado = data_in + caldays(0:(length(t_dado)-1));

Nw_i=diff(Nwm);
Nw_d=diff(Dm);

r_nwi = (Nw_cases)-Nw_i(1:length(Nw_cases))';
r_nwd = (Nw_deaths)-Nw_d(1:length(Nw_deaths))';
r_h = H_dado-Hm(1:length(H_dado))';
r_u = U_dado-Um(1:length(U_dado))';

SSE_nwi = norm(r_nwi).^2;
SSE_nwd = norm(r_nwd).^2;
SSE_h = norm(r_h).^2;
SSE_u = norm(r_u).^2;

SST_nwi = norm(Nw_cases-mean(Nw_i(1:length(Nw_cases))))^2;
SST_nwd = norm(Nw_deaths-mean(Nw_d(1:length(Nw_deaths))))^2;
SST_h = norm(H_dado-mean(Hm(1:length(H_dado))))^2;
SST_u = norm(U_dado-mean(Um(1:length(U_dado))))^2;

% R2 - correlation factor
R2_nwi = 1 - SSE_nwi/SST_nwi;
R2_nwd = 1 - SSE_nwd/SST_nwd;
R2_h = 1 - SSE_h/SST_h;
R2_u = 1 - SSE_u/SST_u;

% Error in %
r_nwi = (abs(r_nwi)./(Nw_cases))*100;
r_nwd = (abs(r_nwd)./Nw_deaths)*100;
r_h = (abs(r_h)./H_dado)*100;
r_u = (abs(r_u)./U_dado)*100;

%% ============ Plot ============

figure(1)
subplot(2,2,3)
H_min = Hm*0.95;
H_max = Hm*1.05;
shadedplot(t_p,H_min',H_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p,Hm,'b','LineWidth',3);
p4=plot(t_dado,H_dado,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(tsim)])
ylabel('Clinical beds occupancy')
      
figure(1)
subplot(2,2,4)
U_min = Um*0.95;
U_max = Um*1.05;
shadedplot(t_p,U_min',U_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p,Um,'b','LineWidth',3);
p4=plot(t_dado,U_dado,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(tsim)])
ylabel('ICU beds occupancy')
    
% figure(3)
% Ic_min = (Nwm)*0.8;
% Ic_max = (Nwm)*1.2;
% shadedplot(t_p,Ic_min',Ic_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
% hold on
% p3=plot(t_p,Nwm,'b','LineWidth',3);
% p4=plot(t_dado,Cases_dado,'k*','LineWidth',2);
% legend([p3 p4],{'SEIIHURD+\Psi Model','Official Data - SESAB'})
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% xlim([t_dado(1) t_dado(tsim)])
% xlabel('Time [days]')
% ylabel('Cumulative Infected Cases')

% figure(4)
% D_min = Dm*0.8;
% D_max = Dm*1.2;
% shadedplot(t_p,D_min',D_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
% hold on
% p3=plot(t_p,Dm,'b','LineWidth',3);
% p4=plot(t_dado,D_dado,'k*','LineWidth',2);
% legend([p3 p4],{'SEIIHURD+\Psi Model','Official Data - SESAB'})
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% xlim([t_dado(1) t_dado(tsim)])
% xlabel('Time [days]')
% ylabel('Fatal Cases')

figure(1)
subplot(2,2,1)
Nw_i_min = Nw_i*0.95;
Nw_i_max = Nw_i*1.05;
shadedplot(t_p(2:end),Nw_i_min',Nw_i_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p(2:end),Nw_i,'b','LineWidth',3);
p4=plot(t_dado(2:end),Nw_cases,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(2) t_dado(tsim)])
ylabel('New cases')

figure(1)
subplot(2,2,2)
Nw_d_min = Nw_d*0.95;
Nw_d_max = Nw_d*1.05;
shadedplot(t_p(2:end),Nw_d_min',Nw_d_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p(2:end),Nw_d,'b','LineWidth',3);
p4=plot(t_dado(2:end),Nw_deaths,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(2) t_dado(tsim)])
ylabel('New fatal cases')

figure(2)
Vm_min = Vm*0.95;
Vm_max = Vm*1.05;
shadedplot(t_p,Vm_min',Vm_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p,NwVm,'b','LineWidth',3);
p4=plot(t_dado,V_dado,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(tsim)])
ylabel('Vaccinated')

% figure(6)
% R_min = Rm*0.8;
% R_max = Rm*1.2;
% shadedplot(t_p,R_min',R_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
% hold on
% p3=plot(t_p,Rm,'b','LineWidth',3);
% p4=plot(t_dado,R_dado,'k*','LineWidth',2);
% legend([p3 p4],{'SEIIHURD+\Psi Model','Official Data - SESAB'})
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% xlim([t_dado(1) t_dado(tsim)])
% xlabel('Time [days]')
% ylabel('Recovered Cases')

lim10 = 10*ones(length(t_dado),1);
figure(10)
plot(t_dado(2:end),r_nwi,'b',t_dado(2:end),r_nwd,'k+')
hold on
plot(t_dado,r_h,'k-',t_dado,r_u,'b--')
plot(lim10,'--r')
xlabel('dias')
ylabel('Erro Data - Modelo [%]')
xlim([t_dado(2) t_dado(end)])
ylim([0 50])
set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
legend(sprintf('New Cases:R^{2} = %0.4f',R2_nwi),sprintf('New Deaths:R^{2} = %0.4f',R2_nwd),sprintf('H:R^{2} = %0.4f',R2_h),sprintf('U:R^{2} = %0.4f',R2_u)) 
title('SEIIHURD+\Psi and R^2')

%% ============ Save Data ============

% saveas(figure(1),[pwd '/Figures/seiihurd_painel.fig']);
% saveas(figure(10),[pwd '/Figures/output_scen1.fig']);
% save('output_data/seiihurd_psi_data.mat')

%% ============ Model simulation ============
function dydt = covid_odes(t,y,Psif)
N=1;
psi=Psif(floor(t+1));

k = 1/4;
gamma_a = 0.18;
gamma_s = 1/4;
gamma_h = 0.25;
gamma_u = 0.13342706158133355;
mi_u = 0.4;
qsi = 0.53;
h = 0.06287612499693644;
h_v=0.015;
mi_h = 0.15;
ome_h = 0.14;
ome_u = 0.29;
% delta = 0.30906304338495505;
delta = 0.31;
p = 0.2;

% gama_v, Iv, omega_v, miVU, omega_u

if t<20.178
beta=2.1317;    %beta=1.3987731952032998;
beta_v=0;
elseif (t>=28.178-8)&&(t< 72.94)
beta=1.7645;    %0.9614724422279308; 
beta_v=0;
elseif (t>=72.94)&&(t< 148)
beta=1.1281;    %0.6657552424857321; 
beta_v=beta/2;
else 
beta=1;
beta_v=beta/2;
end

%  New parameters
%  Vaccination compartments parameters are 0 until starts the vaccination


if t<100
    tau=0;
    delta_av = 0;
    delta_sv = 0;
    phi_e = 0;
    k_v=0;
    p_v = 0;
    gamma_av = 0;
    gamma_sv = 0;
    gamma_vu=0;
    qsi_v = 0 ;
    eps=0;
    mi_vh=0;
    mi_vu=0;
    gamma_vh=0;
else
    delta_av = 0.1;
    delta_sv = 0.2;
    phi_e = 0.8;
%     tau=(6.6976e-05)*1.1;
    tau=1/100000;
    k_v=1/3;
    p_v = 0.1;
    gamma_av = 1/2;
    gamma_sv = 1/2;
    gamma_vu=0.19;
    qsi_v = 0.99;
    eps=0.7;
    mi_vh=0.07;
    mi_vu=0.25;
    gamma_vh=0.25;
end


S=y(1);
E=y(2);
Sv=y(3);
Ev=y(4);
Ia=y(5);
Is=y(6);
Iav=y(7);
Isv=y(8);
H=y(9);
Hv=y(10);
U=y(11);
Uv=y(12);
R=y(13);
Rv=y(14);
D=y(15);
Nw=y(16);
NwV=y(17);

dSdt = -(1-psi)*beta*S*(Is+delta*Ia+delta_av*Iav+delta_sv*Isv)/N - tau*S;
dEdt = (1-psi)*beta*S*(Is+delta*Ia+delta_av*Iav+delta_sv*Isv)/N - k*E;
dSvdt = tau*S - (1-psi)*beta_v*Sv*(Is+delta*Ia+delta_av*Iav+delta_sv*Isv)/N - phi_e*eps*Sv;
dEvdt = (1-psi)*beta_v*Sv*(Is+delta*Ia+delta_av*Iav+delta_sv*Isv)/N - k_v*Ev;
dIadt = (1-p)*k*E - gamma_a*Ia;
dIsdt = p*k*E - gamma_s*Is;
dIavdt = (1-p_v)*k_v*Ev - gamma_av*Iav;
dIsvdt = p_v*k_v*Ev - gamma_sv*Isv;
dHdt = h*qsi*gamma_s*Is + (1-mi_u+ome_u*mi_u)*gamma_u*U - gamma_h*H;
dHvdt = h_v*qsi_v*gamma_sv*Isv + (1-mi_vu+ome_u*mi_vu)*gamma_vu*Uv - gamma_vu*Hv;
dUdt = h*(1-qsi)*gamma_s*Is + ome_h*gamma_h*H - gamma_u*U;
dUvdt = h_v*(1-qsi_v)*gamma_sv*Isv + ome_h*gamma_h*Hv - gamma_u*Uv;
dRdt = gamma_a*Ia + (1-h)*gamma_s*Is + (1-mi_h)*(1-ome_h)*gamma_h*H - phi_e*eps*Sv;
dRvdt = gamma_av*Iav + (1-h)*gamma_sv*Isv + (1-mi_vh)*(1-ome_h)*(gamma_vh*Hv);
dDdt = (1-ome_h)*(mi_h*gamma_h*H +mi_vh*gamma_vh*Hv) + (1-ome_h)*(mi_u*gamma_u*U+mi_vu*gamma_vu*Uv);
dNwdt = p*k*E + p_v*k_v*Ev;
dNwVdt = tau*S;

dydt = [dSdt; dEdt; dSvdt; dEvdt; dIadt; dIsdt; dIavdt; dIsvdt; dHdt; dHvdt; dUdt; dUvdt; dRdt; dRvdt; dDdt; dNwdt; dNwVdt];

end