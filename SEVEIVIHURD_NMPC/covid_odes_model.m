function dydt = covid_odes_model(t,y,Psif,N0,x)
N=1;
Psi=Psif(floor(t+1));

kappa = 1/4;
gamma_a = 1/3.5;
gamma_s = 1/4;
gamma_h = 0.1776863042741119;
gamma_u = 0.13342706158133355;
mi_u = 0.4;
qsi = 0.53;
h = 0.06287612499693644;
mi_h = 0.15;
ome_h = 0.14;
ome_u = 0.29;
delta = 0.30906304338495505;
p = 0.2;

if t<20.178
beta=2.1317;    %beta=1.3987731952032998;
elseif (t>=28.178-8)&&(t< 72.94)
beta=1.7645;    %0.9614724422279308;    
elseif (t>=72.94)&&(t< 148)
beta=1.1281;    %0.6657552424857321; 
else 
    beta=1;
end

%% Vaccine Parameter

if t<100
    tau=0;
    delta_av = 0;
    delta_sv = 0;
    phi_e = 0;
    k_v=0;
    p_v = 0;
    gamma_av = 0;
    gamma_sv = 0;
    qsi_v = 0 ;
else
    delta_av = 0.30906304338495505/20;
    delta_sv = 0.30906304338495505/20;
    phi_e = 0.8;
    tau=(6.6976e-05)*1.1;
    k_v=0.1;
    p_v = 0.1;
    gamma_av = 1/3.5;
    gamma_sv = 1/4;
    qsi_v = 0.53 ;
end

%% ODES

S=y(1);
E=y(2);
V=y(3);
Ev=y(4);
Ia=y(5);
Is=y(6);
Iav=y(7);
Isv=y(8);
H=y(9);
U=y(10);
R=y(11);
Rv=y(12);
D=y(13);
Nw=y(14);
NwV=y(15);

dSdt = -(1-Psi)*beta*S*(Is+delta*Ia+delta_av*Iav+delta_sv*Isv)/N - tau*S;
dEdt = (1-Psi)*beta*S*(Is+delta*Ia+delta_av*Iav+delta_sv*Isv)/N - kappa*E;
dVdt = tau*S - (1-Psi)*beta*V*(Is+delta*Ia+delta_av*Iav+delta_sv*Isv)/N - phi_e*V;
dEvdt = (1-Psi)*beta*V*(Is+delta*Ia+delta_av*Iav+delta_sv*Isv)/N - k_v*Ev;
dIadt = (1-p)*kappa*E - gamma_a*Ia;
dIsdt = p*kappa*E - gamma_s*Is;
dIavdt = (1-p_v)*k_v*Ev - gamma_av*Iav;
dIsvdt = p_v*k_v*Ev - gamma_sv*Isv;
dHdt = qsi_v*gamma_sv*Isv + h*qsi*gamma_s*Is + (1-mi_u+ome_u*mi_u)*gamma_u*U - gamma_h*H;
dUdt = h*(1-qsi)*gamma_s*Is + ome_h*gamma_h*H - gamma_u*U;
dRdt = gamma_a*Ia + (1-h)*gamma_s*Is + (1-mi_h)*(1-ome_h)*gamma_h*H;
dRvdt = gamma_av*Iav + (1-qsi_v)*gamma_sv*Isv + phi_e*V;
dDdt = (1-ome_h)*mi_h*gamma_h*H + (1-ome_u)*mi_u*gamma_u*U;
dNwdt = p*kappa*E + p_v*k_v*Ev;
dNwVdt = tau*S;

dydt = [dSdt; dEdt; dVdt; dEvdt; dIadt; dIsdt; dIavdt; dIsvdt; dHdt; dUdt; dRdt; dRvdt; dDdt; dNwdt; dNwVdt];

end